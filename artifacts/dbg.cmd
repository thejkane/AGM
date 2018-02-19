launch $b{1} ./performance_test --launcher-args="-b -n 1  -N 1 -d 1" --args="--poll-task 1 --threads 1 --scale 10 --degree 16 --num-sources 3 --coalescing-size 42000 --distribution-coalescing-size 142000 --flush 18 --delta 10 --ds nodeq --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_ds --verify"

launch $b{1} ./performance_test --launcher-args="-b -n 1  -N 1 -d 1" --args="--poll-task 1 --threads 1 --scale 10 --degree 16 --num-sources 3 --coalescing-size 42000 --distribution-coalescing-size 142000 --flush 18 --delta 10 --ds numaq --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_ds --verify"



aprun -b -n 2  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 20 --degree 16 --num-sources 3 --coalescing-size 42000 --distribution-coalescing-size 142000 --flush 18 --delta 10 --ds numaq --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_ds --verify


launch $b{1} ./performance_test --launcher-args="-b -n 1  -N 1 -d 31" --args="--poll-task 1 --threads 31 --scale 15 --degree 16 --num-sources 3 --coalescing-size 42000 --distribution-coalescing-size 142000 --flush 18 --delta 10 --ds numaq --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_ds --verify"

launch $b{2} ./performance_test --launcher-args="-b -n 2  -N 1 -d 31" --args="--poll-task 1 --threads 31 --scale 15 --degree 16 --num-sources 3 --coalescing-size 42000 --distribution-coalescing-size 142000 --flush 18 --delta 10 --ds numaq --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_ds --verify"


aprun -b -n 2  -N 1 -d 31 ./performance_test --poll-task 1 --threads 31 --scale 20 --degree 16 --num-sources 3 --coalescing-size 42000 --distribution-coalescing-size 142000 --flush 18 --delta 10 --ds numaq --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_ds --verify


aprun -b -n 2  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 15 --degree 16 --num-sources 3 --coalescing-size 42000 --distribution-coalescing-size 142000 --flush 18 --ds numaq --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_dc --verify


launch $b{2} ./performance_test --launcher-args="-b -n 2  -N 1 -d 31" --args="--poll-task 1 --threads 31 --scale 15 --degree 16 --num-sources 3 --coalescing-size 42000 --distribution-coalescing-size 142000 --flush 18 --ds numaq --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_dc --verify"
