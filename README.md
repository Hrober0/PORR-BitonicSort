# Bitonic sort

Sekwencyjne sortowanie bitoniczne w C11. Wspólny kod jest w `common/`, program w `seq/`. Długość tablicy to zawsze potęga dwójki: `n = 2^log2_n`.

## Struktura

```
common/          gen_data, sortowanie, scalanie, pomiary, CSV
seq/             wersja sekwencyjna
stage1/omp/      OpenMP
stage1/proc/     procesy (fork + mmap)
stage2/cuda/     CUDA
stage2/opencl/   OpenCL
stage2/mpi/      MPI
stage2/k8s/      Kubernetes
results/         pliki CSV
```

## Uruchomienie

```bash
make
./bitonic_seq <log2_n> <p> <seed> [results.csv] [run_id]
```

Przykład: `./bitonic_seq 20 1 1`. Wersja sekwencyjna zapisuje `p=1`. Domyślny plik wyniku to `results/seq.csv`. `status=OK` oznacza, że norma Cauchy'ego i residual są równe `0`. Czas obejmuje tylko `bitonic_sort_seq`.

Kalibracja rozmiaru, cel około 120 s:

```bash
make calibrate
```

Na maszynie testowej `n_max` to `log2_n = 27` (około 84 s). `2^28` zajmuje 165–190 s.

## Kolejne kroki

Wersje równoległe i rozproszone korzystają z `common/`: blok sortuje `bitonic_sort_seq`, scalenie robi `bitonic_merge_seq`, dane pochodzą z `gen_data` (ten sam seed daje tę samą tablicę), a czas, sprawdzenie wyniku i zapis CSV biorą się z gotowych funkcji. OpenMP i procesy różnią się tylko uruchomieniem workerów i synchronizacją.
