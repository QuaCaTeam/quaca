# Parallelization {docsify-ignore-all}

## Compute friction in parallel

Using the default main function of *QuaCa* a simple parallelization of the for-loop was implemented using the [OpenMP API](https://www.openmp.org/resources/refguides/). By default all avaiable threads are used. To limit the number of threads a flag can be set
``` bash
quaca/bin/./Friction --file ../data/MyInputFile.json --threads number_of_threads
```
where `number_of_threads` represents the number of threads that you want to use.

## Parallelize your own function

### Parallelized for loop

In most calcluations using  *QuaCa* we will implement a for-loop iterating through some quantity in the formula. In these cases a parallelization can be implemented in straight forward fashion. First of all we have to start a parallel region by the macro
``` cpp
#pragma omp parallel
{
  /* Your code */
}
```
Every variable defined before this macro is shared between all the threads. For every variable defined within the brackets each threads creates their own copy. Therefore, any variable for which variables/values need to be changed indepedent of the other threads have to be defined after the macro.
Another possibility is to define the variable to be private
``` cpp
int x = 0;
int y = 0;
#pragma omp parallel private(y)
{
 int z = 0;
 std::cout << x++ << "\t" << y++ << "\t" << z++ << std::endl;
} 
```
In this example the variable `int x` gets increased by each thread independently resulting in a total value of `x=2`. The variable y is declared private, therefore each thread creates it's own copy of y, such that y gets increased just once. Finally the variable z is created within the brackets, such that similar to `int y` each thread has it's own copy of `int z`.
Within this parallel region we can now easily create a parallelized for-loop by another macro
``` cpp
#pragma omp parallel
{
  #pragma omp for
  for(int i = 0; i < some_number; ++i)
  {
    /* Your code */
  }
}
```

### Dynamic scheduling

OpenMp automatically divides all the runs of the for-loop evenly between the different threads. In many cases there are specific values, for which calculations take more time. Therefore, it would be more convenient if each thread would takes a new task, as soon as the last jobs is finished, instead of a static division of all tasks at the beginning of the calculation. This can be achieved by setting `schedule` to be `dynamic`.
``` cpp
#pragma omp parallel
{
  #pragma omp for schedule(dynamic)
  for(int i = 0; i < some_number; ++i)
  {
  /* Your code */
  }
}
```

### Setting the number of threads
So far OpenMP requests all possible threads avaiable on the machine for the computations. In some cases it might be preferable if the number of threads could be set, like in the case of `./Friction` where we can set the number with the flag `--threads some_number`. This can be achieved by setting the variable `num_threads`.
``` cpp
#pragma omp parallel num_threads(number_of_threads)
{
  #pragma omp for schedule(dynamic)
  for(int i = 0; i < some_number; ++i)
  {
  /* Your code */
  }
}
```
Further explanations can be found in the [official documentation](https://www.openmp.org/resources/refguides/0). Another useful tutorial can be found [here](https://helloacm.com/simple-tutorial-with-openmp-how-to-use-parallel-block-in-cc-using-openmp/).
