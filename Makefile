export PATH := $(HOME)/.elan/bin:/opt/homebrew/opt/llvm/bin:/opt/homebrew/bin:/usr/local/bin:$(PATH)

FLOW_DIR ?= ../FLOW
CC       ?= clang
CXX      ?= clang++
LEAN_BIN ?= lean
LEANCHECKER ?= leanchecker

CFLAGS   := -std=c17 -O3 -Wall -Wextra -Werror -Iinclude -I$(FLOW_DIR)/core
RING0_CF := -std=c17 -ffreestanding -nostdlib -fno-builtin -fno-stack-protector \
            -mgeneral-regs-only -fno-jump-tables -Wall -Wextra -Werror -O2 \
            -Iinclude -I$(FLOW_DIR)/core

BUILD    := build
TARGET   := $(BUILD)/robo_cli
RING0_O  := $(BUILD)/robo_ring0.o
OLEAN_DIR := $(BUILD)/lean

FLOW_BIN ?= $(FLOW_DIR)/build/flow
AUDITOR_BIN ?= $(HOME)/.local/bin/auditor

HEADERS := include/robo_dynamics_ring0.h include/robo_cbf_ring0.h include/robo_arm2_ring0.h \
           include/robo_qp_nd_ring0.h include/robo_arm6_ring0.h

.PHONY: all ring0 check-headers formal test demo clean

all: $(TARGET) ring0 check-headers formal

$(BUILD):
	mkdir -p $(BUILD)

$(OLEAN_DIR):
	mkdir -p $(OLEAN_DIR)

$(RING0_O): src/robo_ring0.c $(HEADERS) $(FLOW_DIR)/core/flow.h $(FLOW_DIR)/core/flow_opt_params.h $(FLOW_DIR)/core/flow_cholesky_sdd.h | $(BUILD)
	$(CC) $(RING0_CF) -c src/robo_ring0.c -o $@

ring0: $(RING0_O)
	@echo "=== [ROBO Ring-0 Gate] Pure-GPR Freestanding Audit (0 undef, 0 udiv, 0 FPU/SIMD, 0 br/blr) ==="
	@test -z "$$(nm -u $(RING0_O))" || (echo "ERROR: Undefined external symbols in $(RING0_O):" && nm -u $(RING0_O) && exit 1)
	@! objdump -d $(RING0_O) | grep -Eq '[[:space:]](udiv|sdiv|divq?|idivq?|blr|br)[[:space:]]' || \
		(echo "ERROR: Division or indirect branch/call leaked into $(RING0_O)!" && objdump -d $(RING0_O) && exit 1)
	@! objdump -d $(RING0_O) | grep -Eq '([[:space:]]|%)(fadd|fsub|fmul|fdiv|scvtf|ucvtf|fcvtzs|fcvtzu|fcmp|fmov|addsd|subsd|mulsd|divsd|cvtsi2sd|cvttsd2si|xmm[0-9]|ymm[0-9]|zmm[0-9])' || \
		(echo "ERROR: FPU/SIMD instruction leaked into $(RING0_O)!" && objdump -d $(RING0_O) && exit 1)
	@echo "  [PASS] $(RING0_O): zero undefined symbols, zero udiv/sdiv/blr, zero FPU/SIMD"

check-headers:
	$(CC) -std=c17 -Wall -Wextra -Werror -Wno-unused-function -I$(FLOW_DIR)/core -Iinclude -fsyntax-only $(HEADERS)
	$(CXX) -std=c++20 -Wall -Wextra -Werror -Wno-unused-function -iquote $(FLOW_DIR)/core -iquote include -x c++-header -fsyntax-only $(HEADERS)
	@echo "  [PASS] All 5 Ring-0 headers clean under C17 + C++20"

$(TARGET): src/robo_cli.c $(HEADERS) $(FLOW_DIR)/core/flow.h $(FLOW_DIR)/core/flow_cholesky_sdd.h | $(BUILD)
	$(CC) $(CFLAGS) src/robo_cli.c -lm -o $@

$(OLEAN_DIR)/RoboSafety.olean: lean/RoboSafety.lean | $(OLEAN_DIR)
	$(LEAN_BIN) -DwarningAsError=true --root=lean -o $@ $<

formal: $(RING0_O)
	@echo "=== [ROBO FLOW v10 + Lean 4 Gate] ArenaSchema, Z3 Falsification, Lean 4 & NSW-0 ==="
	@mkdir -p $(BUILD) $(OLEAN_DIR)
	"$(FLOW_BIN)" schema check ir/robo_wire.fspec
	"$(FLOW_BIN)" schema emit c17 ir/robo_wire.fspec $(BUILD)/robo_wire.h
	"$(FLOW_BIN)" schema emit smt2 ir/robo_wire.fspec $(BUILD)/robo_wire_schema.smt2
	@test "$$(z3 -T:5 $(BUILD)/robo_wire_schema.smt2 | head -n 1)" = "unsat"
	@z3 -T:5 proofs/robo_falsify.smt2 > $(BUILD)/robo_falsify.out
	@! grep -qiE "(timeout|unknown|unsupported|\(error)" $(BUILD)/robo_falsify.out
	@test "$$(grep -cxE "sat|unsat" $(BUILD)/robo_falsify.out)" -eq 8
	@if command -v $(LEAN_BIN) >/dev/null 2>&1 && command -v $(LEANCHECKER) >/dev/null 2>&1; then \
		$(MAKE) $(OLEAN_DIR)/RoboSafety.olean && \
		LEAN_PATH=$(OLEAN_DIR) $(LEANCHECKER) RoboSafety && \
		! grep -Eq '\bsorry\b|\bnative_decide\b|\baxiom\b' lean/RoboSafety.lean && \
		python3 -c 'import re; txt=open("lean/RoboSafety.lean").read(); thms=re.findall(r"^\s*theorem\s+(\S+)", txt, re.M); assert len(thms)>=8, f"Too few theorems: {len(thms)}"; open("$(BUILD)/AxiomAudit.lean","w").write("import RoboSafety\n" + "".join(f"#print axioms RoboSafety.{t}\n" for t in thms)); print(f"  [LEAN4] Auditing {len(thms)} theorems...")' && \
		ax_out="$$(LEAN_PATH=$(OLEAN_DIR) $(LEAN_BIN) -DwarningAsError=true $(BUILD)/AxiomAudit.lean)" && \
		! printf '%s\n' "$$ax_out" | grep -Eq 'sorryAx|Lean\.ofReduceBool|Lean\.trustCompiler' && \
		echo "  [PASS] lean/RoboSafety.lean compiled + leanchecker verified (0 sorryAx, 0 native_decide)"; \
	else \
		echo "  [SKIP] lean/leanchecker not found on PATH"; \
	fi
	@if [ -x "$(AUDITOR_BIN)" ]; then "$(AUDITOR_BIN)" bug src/robo_ring0.c -Iinclude -I$(FLOW_DIR)/core --check-nsw --strict; fi

test: $(TARGET) ring0 check-headers formal
	./$(TARGET) selftest
	python3 oracles/robo_rk4_oracle.py --cli ./$(TARGET) --episodes 2

demo: test
	./$(TARGET) demo --report ROBO_VERIFICATION_REPORT.md

clean:
	rm -rf $(BUILD)
