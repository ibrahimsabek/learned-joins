# learned-joins

This is the source code and running commands that are used for "Revisiting In-Memory Join Algorithms with Learned Models" submission. A configuration file is provided to build this project with CMake. In the project directory, run:

```
mkdir -p build/release
cd build/release
```

All configurations of join algorithms and generating/loading the synthetic datasets should be set using two "sh" scripts: `scripts/evaluation/base_configs_maker.sh` and `scripts/evaluation/eth_configs_maker.sh`. These scripts will also generate two configuration files: `include/configs/base_configs.h` and `include/configs/eth_configs.sh`. Examples of setting the parameters in these scripts exist in the scripts we used to run our experiments `scripts/evaluation/*.sh`

You can run any of the existing experiment scripts, after setting the approperiate configurations, from the `release` folder using the following command:

```
sh ../../scripts/evaluation/learned_imv_sort_join_runner.sh 
```


To download the [SOSD benchmark datasets](https://github.com/learnedsystems/SOSD), use the script file `scripts/sosd/sosd.sh`. 


