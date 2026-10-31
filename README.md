# Collatz odd cores of Mersenne numbers: reproducibility code

Code accompanying

> W. Zhou, *A 2-adic logarithm law and analytic isometries for the Collatz odd cores of Mersenne numbers*.

The script `reproduce.py` recomputes every table and every numerical statement
in the paper and compares each value with the value printed in the paper.

## Requirements
Python ≥ 3.8, standard library only. All computations use exact integer arithmetic,
except a few floating-point approximations of constants that the paper also gives
only approximately (h₂(α), ρ, log₂ ρ, log₂(α/(1−α))).

## Usage
```bash
python reproduce.py                 # CSV files go to ./output
python reproduce.py --outdir results
```
The script prints `PASS`/`FAIL` for each check and exits with status 0 if all checks pass.
A full report is written to `output/check_report.csv`.

## What is reproduced

| Paper | Output file | Function |
|---|---|---|
| Table 2 (b_L = λ mod 2^L), λ ≡ 573 mod 1024, 48 bits of λ | `table2_lambda_mod_2L.csv`, `lambda_binary_digits.csv` | `section_2adic` |
| Theorems 3.3, 3.4, Cor. 3.5, 3.6, 3.9(iii), Prop. 3.7 | (checks) | `section_2adic` |
| Table 3 (Terras exceptional sets), Theorems 4.4, 4.5(i), integer inequalities | `table3_exceptional_sets.csv` | `section_terras` |
| Table 4, Theorem 5.1 (67 crossing n ≤ 800, density 21/128), Remark 5.2, Thm 5.4, Cor 5.5 | `table4_one_step_crossing.csv`, `crossing_n_6_to_800.csv` | `section_crossing` |
| Tables 5, 6, A2, A3, A4; Theorems 6.2, 6.5, 6.7, 6.11, 6.15; Cor. 6.12, 6.13, 6.16 | `table5_choi_rk.csv`, `table6_beta_v_k4.csv`, `tableA2_block_parameters.csv`, `tableA3_Fj_factorization.csv`, `tableA4_pratt_certificates.csv` | `section_ternary` |
| Table A1, Theorems 2.2, 2.4, 2.5, Lemma 2.1, 2.3, Prop. 2.6, 2.7 | `tableA1_tau_Mn.csv` (2 ≤ n ≤ 100) | `section_prefix` |
| Table A5 and the statistics of first meetings (149 n, max 1472 steps, target frequencies) | `tableA5_first_meetings.csv`, `meetings_all_even_n_4_to_300.csv` | `section_meetings` |

### Remarks on the methods
* λ = Log(9)/8 is computed modulo 2^117 from the series Log(1+8) = Σ (−1)^{j+1} 8^j / j.
* B_q mod 2^L and Φ_q(r) mod 2^P are computed by modular exponentiation of 3, without forming 3^{2^q} in full.
* |E_L^T| is computed by dynamic programming over parity vectors and cross-checked by brute force for L ≤ 14.
* F_4 and F_5 are factored by trial division over candidates ≡ 1 (mod 2^{j+1}) (Theorem 3.9 of the paper), and the Pratt certificates are checked recursively.
* Theorem 6.11 is checked digit by digit: the ternary string of T₁^{k+v₂(y)}(y·3^N) is compared with head + copies of β_v + prefix + tail.

## License
MIT (see `LICENSE`).
# mersennecollatz2adic
mersennecollatz2adic code
