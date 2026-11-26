export PATH := $(HOME)/.elan/bin:/opt/homebrew/opt/llvm/bin:/opt/homebrew/bin:/usr/local/bin:$(PATH)

FLOW_DIR ?= ../FLOW
CC       ?= clang
CXX      ?= clang++

CFLAGS   := -std=c17 -O3 -Wall -Wextra -Werror -Iinclude -I$(FLOW_DIR)/core
RING0_CF := -std=c17 -ffreestanding -nostdlib -fno-builtin -fno-stack-protector \
            -mgeneral-regs-only -fno-jump-tables -Wall -Wextra -Werror -O2 \
            -Iinclude -I$(FLOW_DIR)/core

BUILD    := build
TARGET   := $(BUILD)/robo_cli
RING0_O  := $(BUILD)/robo_ring0.o

FLOW_BIN ?= $(FLOW_DIR)/build/flow
AUDITOR_BIN ?= $(HOME)/.local/bin/auditor

.PHONY: all ring0 check-headers formal test demo clean

all: $(TARGET) ring0 check-headers formal

$(BUILD):
	mkdir -p $(BUILD)

$(RING0_O): src/robo_ring0.c include/robo_dynamics_ring0.h include/robo_cbf_ring0.h include/robo_arm2_ring0.h $(FLOW_DIR)/core/flow.h $(FLOW_DIR)/core/flow_opt_params.h | $(BUILD)
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
	$(CC) -std=c17 -Wall -Wextra -Werror -Wno-unused-function -I$(FLOW_DIR)/core -Iinclude -fsyntax-only include/robo_dynamics_ring0.h include/robo_cbf_ring0.h include/robo_arm2_ring0.h
	$(CXX) -std=c++20 -Wall -Wextra -Werror -Wno-unused-function -iquote $(FLOW_DIR)/core -iquote include -x c++-header -fsyntax-only include/robo_dynamics_ring0.h include/robo_cbf_ring0.h include/robo_arm2_ring0.h
	@echo "  [PASS] robo_dynamics_ring0.h + robo_cbf_ring0.h + robo_arm2_ring0.h clean under C17 + C++20"

$(TARGET): src/robo_cli.c include/robo_dynamics_ring0.h include/robo_cbf_ring0.h include/robo_arm2_ring0.h $(FLOW_DIR)/core/flow.h | $(BUILD)
	$(CC) $(CFLAGS) src/robo_cli.c -lm -o $@

formal: $(RING0_O)
	@echo "=== [ROBO FLOW v10 Gate] ArenaSchema, Z3 ECBF Falsification & NSW-0 ==="
	@mkdir -p $(BUILD)
	"$(FLOW_BIN)" schema check ir/robo_wire.fspec
	"$(FLOW_BIN)" schema emit c17 ir/robo_wire.fspec $(BUILD)/robo_wire.h
	"$(FLOW_BIN)" schema emit smt2 ir/robo_wire.fspec $(BUILD)/robo_wire_schema.smt2
	@test "$$(z3 -T:5 $(BUILD)/robo_wire_schema.smt2 | head -n 1)" = "unsat"
	@z3 -T:5 proofs/robo_falsify.smt2 > $(BUILD)/robo_falsify.out
	@! grep -qiE "(timeout|unknown|unsupported|\(error)" $(BUILD)/robo_falsify.out
	@test "$$(grep -cxE "sat|unsat" $(BUILD)/robo_falsify.out)" -eq 4
	@if [ -x "$(AUDITOR_BIN)" ]; then "$(AUDITOR_BIN)" bug src/robo_ring0.c -Iinclude -I$(FLOW_DIR)/core --check-nsw --strict; fi

test: $(TARGET) ring0 check-headers formal
	./$(TARGET) selftest

demo: test
	./$(TARGET) demo --report ROBO_VERIFICATION_REPORT.md

clean:
	rm -rf $(BUILD)
