# Hot loop comparison: BLAST+ `BlastAaWordFinder` and the Java C2 `WordFinder`

Observations only. Addresses are as printed in the listings named below, so each line can be checked against its source.

## Sources

| Side | Code | Listing | Notes |
|---|---|---|---|
| BLAST+ | `blastp` 2.17.0 (`~/opt/ncbi-blast-2.17.0+/bin/blastp`, stripped) | `scratchpad/blastp-wf-annot.txt`: `perf annotate` of `BlastAaWordFinder`, 62,575 samples, `instructions:u` | The query and database of that `perf` run are not recorded in the file. Compiler and flags of the release binary were not checked. |
| Java | `WordFinder.find` on `main` (standalone compile, 169 bytes), JDK 25.0.4 (Temurin), C2 | `scratchpad/find-now.s` (from `asm-bench.log`, the listing of the compiled method, `WordFinder::find`) | No sample counts for this listing. |
| Java | `WordFinder.extendLeft` (unchanged on `main` since commit `ca51503`) | `scratchpad/asm-exp6.log`, addresses `0x7f703bef7358`–`0x7f703bef7397` | Same JDK. No sample counts. |

Two regions were compared: the neighbour loop with the diagonal update (loop 1), and the X-drop extension loop (loop 2). "Usual hit" below means the path taken when the diagonal entry is not flagged and the last hit is a window or more back.

---

## Loop 1: neighbour scan and diagonal update

### BLAST+ (`blastp-wf-annot.txt`)

Usual hit: `ff3090` → `ff30cb` → `ff3030` → `ff3040` → `ff3170` → `ff3182` → back to `ff3090`. 30 instructions, 5 with stack operands, 0 stack stores to the frame.

| Address | Instruction | Sample % | Observation |
|---|---|---|---|
| ff3090 | `movq (%rsp),%rax` | 3.73 | Loads the loop end from the frame (`(%rsp)`). |
| ff3094 | `addq $0x8,%r13` | 0.29 | Advances `%r13` by 8: each neighbour entry is 8 bytes. |
| ff3098 | `cmpq %rax,%r13` | 0.82 | Loop bound test. |
| ff309b | `je ff30f2` | 0.00 | Exit, not taken. |
| ff309d | `movl (%r13),%r12d` | 0.50 | First int of the entry into `%r12d`. |
| ff30a1 | `movl 0x4(%r13),%ebx` | 1.03 | Second int of the entry into `%ebx`. |
| ff30a5 | `movl 0x14(%rsp),%eax` | 0.15 | Reloads a value from the frame into `%eax` (the diagonal mask, per the `andl` below). |
| ff30a9 | `movq 0x8(%rsp),%rsi` | 0.95 | Reloads the diagonal array base from the frame into `%rsi`. |
| ff30ae | `movl %r12d,%r15d` | 0.42 | Copies `%r12d` to `%r15d`. |
| ff30b1 | `subl %ebx,%r15d` | 0.83 | `%r15d` = first int − second int. |
| ff30b4 | `andl %r15d,%eax` | 0.20 | `%eax` = `%eax` & `%r15d`. |
| ff30b7 | `cltq` | 1.15 | Sign-extends `%eax` to 64 bits. |
| ff30b9 | `leaq (%rsi,%rax,4),%rbp` | 0.55 | `%rbp` = address of the diagonal entry (4-byte entries). |
| ff30bd | `movl (%rbp),%eax` | 1.38 | Loads the diagonal entry into `%eax`. |
| ff30c0 | `leal (%rax,%rax),%r14d` | **13.54** | `%r14d` = 2 × entry. The largest sample count in the function. The sample percentage is on this instruction, which runs right after the load above. |
| ff30c4 | `sarl %r14d` | 2.17 | Arithmetic shift right of `%r14d` by 1: `%r14d` = entry with bit 31 copied into bit 30. |
| ff30c7 | `cmpb $0x0,0x3(%rbp)` | 2.14 | Tests the sign of the entry (its top byte). The flag is the sign bit. |
| ff30cb | `jns ff3030` | 7.84 | Not flagged: go to `ff3030`. |
| ff3030 | `movl 0x10(%rsp),%esi` | 1.35 | Reloads a value from the frame into `%esi`. |
| ff3034 | `movl %ebx,%edx` | 0.16 | `%edx` = `%ebx`. |
| ff3036 | `subl %esi,%r14d` | 0.44 | `%r14d` = `%r14d` − `%esi`. |
| ff3039 | `subl %r14d,%edx` | 0.41 | `%edx` = `%edx` − `%r14d`. |
| ff303c | `cmpl %edx,0x18(%rsp)` | 1.91 | Compares `%edx` with a frame value. |
| ff3040 | `jle ff3170` | 1.19 | Taken on the usual hit (the path continues at `ff3170`). |
| ff3170 | `addl %esi,%ebx` | 1.63 | `%ebx` = `%ebx` + `%esi`. |
| ff3172 | `andl $0x80000000,%eax` | 0.81 | Keeps only the flag bit of the old entry. |
| ff3177 | `andl $0x7fffffff,%ebx` | 0.83 | Clears the flag bit of the new value. |
| ff317d | `orl %ebx,%eax` | 0.07 | Combines the flag and the new value. |
| ff317f | `movl %eax,(%rbp)` | 0.78 | Stores the new entry. |
| ff3182 | `jmp ff3090` | 0.49 | Back to the loop head. |

The non-zero samples that are not on the usual path, for reference: `ff3046` (5.43%, `cmpl %edx,0x20(%rsp)`), `ff3090`, and the rare-path call to `BSearchContextInfo` at `ff3057`.

### Java C2 (`find-now.s`)

Usual hit, instruction addresses as printed: `0x9a20`–`0x9a8d`, `0x9a93`–`0x9aad`, `0x9b7e`–`0x9ba3`. 38 instructions, 16 with stack operands, 4 other memory operands.

| Address | Instruction | Observation |
|---|---|---|
| 0x9a20 | `movl 0x10(%rax,%r10,4),%r11d` | Loads `offsets[k]` into `%r11d`. |
| 0x9a25 | `movl %r11d,0x2c(%rsp)` | Stores `%r11d` to the frame. |
| 0x9a2a | `movl %r10d,0x28(%rsp)` | Stores `%r10d` (k) to the frame. |
| 0x9a2f | `movq %rax,8(%rsp)` | Stores `%rax` (offsets array) to the frame. |
| 0x9a34 | `movq 0x68(%rsp),%r10` | Reloads `this` from the frame. |
| 0x9a39 | `movl 0x44(%r10),%r11d` | Loads the `diagState` field into `%r11d`. |
| 0x9a3d | `movl 0x2c(%rsp),%r10d` | Reloads q from the frame. |
| 0x9a42 | `subl 0xa8(%rsp),%r10d` | `%r10d` = q − s (s from the frame). |
| 0x9a4a | `movq 0x68(%rsp),%r8` | Reloads `this`. |
| 0x9a4f | `andl 0x18(%r8),%r10d` | `%r10d` &= `this.diagMask` (field load). |
| 0x9a53 | `movl %r10d,0x30(%rsp)` | Stores the diagonal index to the frame. |
| 0x9a58 | `movl 0xc(%r12,%r11,8),%r10d` | Loads `diagState.length` (bounds check). |
| 0x9a5d | `movl 0x30(%rsp),%r9d` | Reloads the diagonal index. |
| 0x9a62 | `cmpl %r10d,%r9d` | Bounds check. |
| 0x9a65 | `jae 0x9d08` | Bounds check, not taken. |
| 0x9a6b | `leaq (%r12,%r11,8),%rcx` | Address of the `diagState` array data. |
| 0x9a6f | `movl 0x1c(%r8),%r8d` | Loads `this.diagOffset`. |
| 0x9a73 | `movl %r9d,%r10d` | Copies the index. |
| 0x9a76 | `movl 0x10(%rcx,%r10,4),%r9d` | Loads the diagonal state into `%r9d`. |
| 0x9a7b | `movl %r8d,%r10d` | Copies `diagOffset`. |
| 0x9a7e | `addl 0xa8(%rsp),%r10d` | `%r10d` = s + diagOffset (s from the frame). |
| 0x9a86 | `testl $0x40000000,%r9d` | Tests the EXTENDED bit (bit 30). |
| 0x9a8d | `jne 0x9ba8` | Not taken on the usual hit. |
| 0x9a93 | `subl %r8d,%r9d` | `%r9d` = state − diagOffset (last). |
| 0x9a96 | `movq 0x68(%rsp),%r11` | Reloads `this`. |
| 0x9a9b | `movl 0xc(%r11),%r8d` | Loads `this.window`. |
| 0x9a9f | `movl 0xa8(%rsp),%r11d` | Reloads s. |
| 0x9aa7 | `subl %r9d,%r11d` | `%r11d` = s − last (diff). |
| 0x9aaa | `cmpl %r8d,%r11d` | diff vs window. |
| 0x9aad | `jge 0x9b7e` | Taken on the usual hit. |
| 0x9b7e | `movl 0x30(%rsp),%r11d` | Reloads the diagonal index. |
| 0x9b83 | `movl %r10d,0x10(%rcx,%r11,4)` | Stores the new diagonal state. |
| 0x9b88 | `movl 0x28(%rsp),%r10d` | Reloads k. |
| 0x9b8d | `incl %r10d` | k + 1. |
| 0x9b90 | `cmpl 0xac(%rsp),%r10d` | Compares k with `end` (from the frame). |
| 0x9b98 | `jge 0x990d` | Exit of the neighbour loop when k reaches `end`. |
| 0x9b9e | `movq 8(%rsp),%rax` | Reloads the offsets array. |
| 0x9ba3 | `jmp 0x9a20` | Back to the neighbour loop head. |

Before the usual-hit path, each neighbour goes through this block (no samples available for this listing):

| Address | Observation |
|---|---|
| 0x9955–0x9963 | Loop-head test (`cmpl %r11d,%r10d`; `jge`). |
| 0x9969–0x9972 | Loads the residue at `start + s + 2`. |
| 0x9978–0x997f | Shifts and ORs the word, then `andl $0x7fff`, the word mask (15 bits). |
| 0x9986–0x9989 | Bounds check on `cellStart` (`cmpl %esi,%r9d; jae`). |
| 0x998f, 0x99a8 | Loads `cellStart[word]` and `cellStart[word+1]`. |
| 0x99ad–0x99b0 | `k < end` test. |
| 0x99b2–0x99c4 | Bounds checks on `offsets` (`cmpl %r14d,%r10d; jae`, `decq`, `cmpq %r13,%rbp; jae`). |
| 0x99ca–0x9a0d | **Ten stores to the frame**, one each for `%r8`, `%rbx`, `%r9d`, `%r11d`, `%esi`, `%r14d`, `%r13`, `%rdi`, `%ecx`, `%edx`, to `0x70`, `0x80`, `0x88`, `0x8c`, `0x90`, `0x94`, `0x98`, `0xa0`, `0xa8`, `0xac`(%rsp). These run on every neighbour (the fall-through into `0x9a20`). |
| 0x990d–0x994e | Ten reloads from the frame (`0x70` to `0xa8`), at the top of the outer (subject position) loop, reached after the neighbour loop ends (`jge 0x990d` at `0x9b98`). |

---

## Loop 2: X-drop extension (leftward)

### BLAST+ (`blastp-wf-annot.txt`, `ff3450`–`ff347c`)

Usual step (no new maximum): 13 instructions, 0 stack operands.

| Address | Instruction | Sample % | Observation |
|---|---|---|---|
| ff3450 | `subq $0x1,%rdx` | 3.93 | Decrements the index `%rdx`. |
| ff3454 | `testl %edx,%edx` | 0.14 | |
| ff3456 | `js ff347e` | 0.00 | Exit when the index is negative, not taken. |
| ff3458 | `movzbl (%r10,%rdx),%edi` | 0.22 | Zero-extends a subject residue into `%edi`. |
| ff345d | `movzbl (%rsi,%rdx),%r11d` | 0.69 | Zero-extends a query residue into `%r11d`. |
| ff3462 | `movq (%r12,%rdi,8),%rdi` | 2.74 | Loads a matrix row pointer from `%r12` indexed by the subject residue. |
| ff3466 | `addl (%rdi,%r11,4),%ecx` | 1.37 | `%ecx` (score) += the matrix entry at the row pointer, indexed by the query residue. |
| ff346a | `movl %eax,%edi` | 8.75 | Copies `%eax` (maxScore) to `%edi`. |
| ff346c | `subl %ecx,%edi` | 0.16 | `%edi` = maxScore − score. |
| ff346e | `cmpl %eax,%ecx` | 3.51 | score vs maxScore. |
| ff3470 | `jle ff3479` | 0.00 | Taken when there is no new maximum. |
| ff3479 | `cmpl %edi,%r15d` | 0.73 | Compares `%r15d` (X-drop) with `%edi`. |
| ff347c | `jg ff3450` | 0.00 | Loop back while `%r15d` > `%edi`. |

On a new maximum, three more instructions run (`movl %edx,%r9d`, `movl %ecx,%eax`, `xorl %edi,%edi`), at `ff3472`–`ff3477` (samples 1.92, 0.18, 0.12).

### Java C2 (asm-exp6.log, `0x7f703bef7358`–`0x7f703bef7397`)

Usual step: 16 instructions, 2 with stack operands. No bounds checks in the loop body. No sample counts.

| Address | Instruction | Observation |
|---|---|---|
| 0x7358 | `movl %ebx,%r10d` | Copies `i`. |
| 0x735b | `addl 4(%rsp),%r10d` | `i + sBase`, with sBase loaded from the frame. |
| 0x7360 | `movsbl 0x10(%rdx,%r10),%r10d` | Sign-extends the subject residue (byte array data at +0x10). |
| 0x7366 | `movl %ebx,%r9d` | Copies `i`. |
| 0x7369 | `addl (%rsp),%r9d` | `i + qBase`, with qBase loaded from the frame. |
| 0x736d | `addl 0x10(%rcx,%r9,4),%r10d` | Adds `rowOfQuery[i + qBase]` (int array, data at +0x10). |
| 0x7372 | `andl $0x3ff,%r10d` | Masks the index to 1023. |
| 0x7379 | `movsbl 0x10(%r13,%r10),%r10d` | Sign-extends `PADDED[index]`. `%r13` holds the static table's address. |
| 0x737f | `addl %r10d,%edi` | score += matrix entry. `%edi` is score. |
| 0x7382 | `cmpl %eax,%edi` | score vs maxScore (`%eax`). |
| 0x7384 | `jg 0x7483` | New maximum: branch out of the loop body to `0x7483`. |
| 0x738a | `cmpl %r8d,%edi` | score vs floor (`%r8d`). |
| 0x738d | `jle 0x7464` | Exit when score ≤ floor. |
| 0x7393 | `decl %ebx` | i − 1. |
| 0x7395 | `cmpl %ebp,%ebx` | i vs the loop bound (`%ebp`). |
| 0x7397 | `jg 0x7358` | Loop back. |

Register roles as used in this loop: `%ebx` = i, `%edi` = score, `%eax` = maxScore, `%r8d` = floor, `%ebp` = loop bound, `%rdx` = subject residue array, `%rcx` = `rowOfQuery` array, `%r13` = static `PADDED` array.

---

## Register use in the hot paths

### Loop 1, BLAST+ (usual hit)

| Register | Use in the listing |
|---|---|
| `%r13` | Neighbour entry pointer (+8 per entry). |
| `%r12d`, `%ebx` | The two ints of the entry (`(%r13)`, `4(%r13)`). |
| `%r15d` | Set to `%r12d` − `%ebx`, then used as the diagonal index before masking. |
| `%eax` | Mask, then diagonal state; then the flag bit. |
| `%rsi` | Diagonal array base, reloaded from `8(%rsp)`. |
| `%rbp` | Address of the diagonal entry. |
| `%r14d` | Entry × 2, then the shifted value, then the difference. |
| `%esi`, `%edx` | Reloaded from the frame (`0x10(%rsp)`), then the difference. |
| `%rax` (loop end) | Reloaded from `(%rsp)` at each iteration. |

Stack operands on the usual path: `(%rsp)`, `0x14(%rsp)`, `0x8(%rsp)`, `0x10(%rsp)`, `0x18(%rsp)`. No stores to the frame.

### Loop 1, Java C2 (usual hit)

| Register | Use in the listing |
|---|---|
| `%r11d` | offsets[k] (q), then `diagState` field, then the diagonal state store index. |
| `%r10d` | Temporary: k, q − s, s + diagOffset, the store value. |
| `%r9d` | Diagonal index, then the state. |
| `%r8` / `%r8d` | `this`, then `diagOffset`, then `window`. |
| `%rcx` | Base of the `diagState` array data. |
| `%r12` | Heap base for compressed references: array and object references are addressed as `(%r12,reg,8)`, e.g. `diagState` at `0xc(%r12,%r11,8)`. |
| `%rax` | Base of the offsets array, spilled to `8(%rsp)` on each neighbour. |
| `%rbp` | Not used on the usual path. |
| Frame | `this` (0x68), q (0x2c), k (0x28), diagonal index (0x30), s (0xa8), end (0xac), offsets array (0x8, and reloaded at 0x9b9e). |

Stack operands on the usual path: 16 (see the table above). Stores to the frame on the usual path: `0x2c`, `0x28`, `8`, `0x30` (4 stores, plus the 10 stores of the neighbour block above).

### Loop 2, BLAST+ (extension step)

| Register | Use in the listing |
|---|---|
| `%rdx` | Index, decremented each step. |
| `%r10` | Subject residue base. |
| `%rsi` | Query residue base. |
| `%r12` | Matrix row-pointer table. |
| `%rdi` | Residue, then the row pointer. |
| `%r11d` | Query residue. |
| `%ecx` | score. |
| `%eax` | maxScore. |
| `%edi` | maxScore − score, then the new-maximum flag. |
| `%r15d` | X-drop (compared with `%edi`). |
| `%r9d` | Best index (written on a new maximum). |

No stack operands in the loop.

### Loop 2, Java C2 (extension step)

| Register | Use in the listing |
|---|---|
| `%ebx` | i. |
| `%edi` | score. |
| `%eax` | maxScore. |
| `%r8d` | floor. |
| `%ebp` | loop bound. |
| `%rdx` | Subject residue array. |
| `%rcx` | `rowOfQuery` array. |
| `%r13` | `PADDED` array. |
| `%r10`, `%r9` | Temporaries (indices and loaded values). |
| `4(%rsp)` | sBase. |
| `(%rsp)` | qBase. |

Two stack operands per step.

---

## Compilation differences

1. **Loop shape.** BLAST counts an index down and tests it with `js`, through a `subq`. Java counts `i` down with `decl` and compares against a register bound (`cmpl %ebp,%ebx`; `jg`).
2. **Matrix lookup.** BLAST loads a row pointer from an 8-byte table indexed by the subject residue (`movq (%r12,%rdi,8)`), then adds the int entry at `(%rdi,%r11,4)`, indexed by the query residue. Java adds a precomputed per-position `rowOfQuery` value to the subject residue, masks with `0x3ff`, and loads from the static `PADDED` table (`movsbl`).
3. **Residue extension.** BLAST zero-extends residues (`movzbl`). Java sign-extends them (`movsbl`).
4. **Bounds checks.** The Java extension loop has none in its body. Its checks are outside the loop. The Java diagonal path has bounds checks on `diagState` (`cmpl`; `jae`), and the neighbour block has checks on `cellStart` and `offsets`. The BLAST+ listing shows no bounds checks in either loop.
5. **Stack traffic.** The BLAST+ loops read the frame on the usual path (5 operands in loop 1, none in loop 2). The Java extension loop reads two frame slots every step (`4(%rsp)`, `(%rsp)`). The Java diagonal path has 16 frame operands, plus 10 stores at the start of each neighbour.
6. **Spill stores in the scan.** The BLAST+ scan loop has no stores to the frame. The Java scan loop stores ten values to the frame on each neighbour (`0x99ca`–`0x9a0d`).
7. **Flag encoding.** BLAST+ uses the sign bit of the diagonal entry (`cmpb $0x0,0x3(%rbp)`; `jns`). It keeps that bit with `andl $0x80000000`. Java uses bit 30 (`testl $0x40000000`), and clears it with `andl $0xbfffffff` (and sets it with `orl $0x40000000`).
8. **New maximum in the extension.** BLAST+ handles it with a short branch (`jle ff3479`) over three instructions in the loop body. Java branches out of the loop body (`jg 0x7483`), so the new-maximum code is not in the loop listing.
9. **Neighbour record.** BLAST+ reads two ints per neighbour entry from an 8-byte record (`movl (%r13)`, `movl 0x4(%r13)`, `addq $0x8,%r13`). Java reads one int per neighbour from `offsets[k]`, with a bounds check on `k`.
10. **Frame size.** BLAST+ frame: `subq $0xc8,%rsp` after six pushes. Java frame: `subq $0x40,%rsp` in `extendLeft`, with spill slots up to `0xac(%rsp)` in `find`.
11. **Compare with the X-drop.** BLAST+ computes maxScore − score into a register and compares it with the X-drop register. Java compares the score with a floor register, a value that changes only on a new maximum.

## Caveats

- The Java listings have no sample counts, so the Java side cannot be weighted by time as the BLAST+ side is. Instruction counts for both sides are from the listings, not from `perf`.
- BLAST+ sample percentages come from one `perf` run of an unrecorded query and database.
- The Java listing of loop 1 is the standalone compile of `WordFinder::find`. In the main build `find` is also compiled inline into the benchmark method, which can produce different code.
- The `0x9955` block and the frame reloads are as printed; the reload block at `0x990d` runs after the neighbour loop ends, not on each neighbour.
