# PA4 RMST Testcases

Expected outputs (`.out`) are tracked in git. Input files (`.dat`) are generated
locally and ignored by git to avoid committing large data files.
Run submissions as:

```sh
python3 generate_testcases.py
./RMST testcases/small.dat answer.dat
diff -u testcases/small.out answer.dat
```

## Files

| Case | Points | Purpose | Expected weight |
| --- | ---: | --- | ---: |
| `small` | 5 | Assignment sample | 13 |
| `mid` | 40,004 | 200 x 200 grid plus duplicate corner points; includes negative coordinates | 280789 |
| `large` | 1,000,000 | Larger grid for performance and memory behavior | 1098900000 |
| `max` | 10,000,000 | Maximum `n`; expected weight requires signed 64-bit integer storage | 999999900000 |

## Grid Answer Formula

For a `cols x rows` grid with horizontal spacing `dx` and vertical spacing `dy`,
the RMST weight is:

```text
if dx <= dy:
    rows * (cols - 1) * dx + (rows - 1) * dy
else:
    cols * (rows - 1) * dy + (cols - 1) * dx
```

Duplicate points add zero-weight edges and do not change the expected RMST weight.

Generate or regenerate all `.dat` files with:

```sh
python3 generate_testcases.py
```
