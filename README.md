# Collatz odd cores of Mersenne numbers: reproducibility code

Code accompanying

> W. Zhou, *A 2-adic logarithm law and analytic isometries for the Collatz odd cores of Mersenne numbers*.

The script `reproduce.py` recomputes every table and every numerical statement in the paper and compares each value with the value printed in the paper.

## Requirements

Python >= 3.8, standard library only. All computations use exact integer arithmetic, except a few floating-point approximations of constants that the paper also gives only approximately (h2(alpha), rho, log2 rho, log2(alpha/(1-alpha))).

## Usage

```
python reproduce.py                 # CSV files go to ./output
python reproduce.py --outdir results
```

The script prints `PASS`/`FAIL` for each check and exits with status 0 if all checks pass.

**The tables of the paper are not printed to the console.** They are written as
machine-readable CSV files into the output directory (created automatically if it
does not exist; default `./output`). Each CSV has a header row whose columns
mirror the columns of the corresponding table in the paper. A reviewer should
open the output directory in Excel or a text editor to inspect the tables.

A successful run ends with a line such as

```
Summary: 75 PASS, 0 FAIL, total time 0.7 s
```

A full report of all checks is written to `output/check_report.csv`.

## Where each table of the paper is generated

| Paper table | Generated CSV file | Function in `reproduce.py` |
|---|---|---|
| Table 2 (b_L = lambda mod 2^L) | `output/table2_lambda_mod_2L.csv` | `section_2adic` |
| 48 binary digits of lambda | `output/lambda_binary_digits.csv` | `section_2adic` |
| Table 3 (Terras exceptional sets) | `output/table3_exceptional_sets.csv` | `section_terras` |
| Table 4 (one-step crossing) | `output/table4_one_step_crossing.csv` | `section_crossing` |
| crossing n, 6 <= n <= 800 | `output/crossing_n_6_to_800.csv` | `section_crossing` |
| Table 5 (Choi's strings r_k) | `output/table5_choi_rk.csv` | `section_ternary` |
| Table 6 (beta_v, k = 4) | `output/table6_beta_v_k4.csv` | `section_ternary` |
| Table A1 (tau(M_n), 2 <= n <= 100) | `output/tableA1_tau_Mn.csv` | `section_prefix` |
| Table A2 (block parameters of B_q) | `output/tableA2_block_parameters.csv` | `section_ternary` |
| Table A3 (factorizations of F_j) | `output/tableA3_Fj_factorization.csv` | `section_ternary` |
| Table A4 (Pratt certificates) | `output/tableA4_pratt_certificates.csv` | `section_ternary` |
| Table A5 (first meetings with earlier cores) | `output/tableA5_first_meetings.csv` | `section_meetings` |
| all meetings, even 4 <= n <= 300 | `output/meetings_all_even_n_4_to_300.csv` | `section_meetings` |
| report of all checks | `output/check_report.csv` | `main` |

All remaining checks (Theorems 2.2, 2.4, 2.5, 3.3, 3.4, 4.4, 4.5, 5.1, 5.4, 6.2,
6.5, 6.7, 6.11, 6.15; Corollaries 3.5, 3.6, 3.9(iii), 5.5, 6.12, 6.13, 6.16;
Propositions 2.6, 2.7, 3.7; Lemmas 2.1, 2.3, 2.9, 6.6; Remarks 5.2, 6.8) are
congruence / inequality / identity verifications; they do not produce a table and
are reported only on the console and in `check_report.csv`.

## How the checks work

For every table, `reproduce.py` recomputes the entries from scratch and compares
them with the values printed in the paper (the `EXP_*` constants at the top of
the script). A `PASS` therefore means that the recomputation reproduces the
published numbers; the CSV files in the output directory are the recomputed
tables themselves. The remaining checks verify the stated congruences, equalities
and inequalities directly over the ranges indicated in the check descriptions.

## Remarks on the methods

- lambda = Log(9)/8 is computed modulo 2^117 from the series Log(1+8) = sum (-1)^{j+1} 8^j / j.
- B_q mod 2^L and Phi_q(r) mod 2^P are computed by modular exponentiation of 3, without forming 3^{2^q} in full.
- |E_L^T| is computed by dynamic programming over parity vectors and cross-checked by brute force for L <= 14.
- F_4 and F_5 are factored by trial division over candidates = 1 (mod 2^{j+1}), and the Pratt certificates are checked recursively.
- Theorem 6.11 is checked digit by digit: the ternary string of T_1^{k+v2(y)}(y * 3^N) is compared with head + copies of beta_v + prefix + tail.

## License

MIT (see `LICENSE`).
