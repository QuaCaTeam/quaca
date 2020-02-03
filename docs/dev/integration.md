# On numerical integration in C++

## The issue
The GSL routine needs to get passed a pointer on a function of the form

```cpp
double my_func(double x, void* p) {
  return x;
}
```

We would like to use a member function of a class, we have defined previously. Somehting like

```cpp
class myClass {
public:
  double my_func(double x, void* p){ return x; };
}
```

Unfortunately the pointer on a globally defined function is not the same, as the pointer to a member functino of a class. This can be understood in the following way. For a general function the compiler knows at compile time (when you do: make) that the function exist and can store the function at some specific location in the memory. This position is constant and we can easily define a pointer pointing to this position. A member function of a class is on the other hand always connected to a specific instance of a class. So this member function only exist after we did something like

```cpp
myClass my_instance();
```

After this instanciation the class has a position in the memory and therefore the member function as well. This happens now at running time and not on compile time. Therefore, both pointers on functions and pointers on member-functions are distinct concepts.

## Possible solutions

### Static member functions

While generally a member function is bound to its instance of the class, a static member function is independent of such an instance, in constrast is has no "this" pointer and therefore cannot access any attribute of the class. Therefore it can also be passed directly to the GSL.

### Void pointer

The void pointer which can be passed in the function for the GSL can generally point to any kind of variable but before it can be used, it has to be casted to the correct type like

```cpp
double my_func(double x, void * p) {
  double* a = static_cast<double*>(p);
  return x+(*a);
}
```

Besides the simple variables like double and int the void pointer can also point to some class or struct. In this class/struct we can now put all variables we need, for example also the pointer to the class, of which we want to integrate a member class from

## The chosen solution

To solve the initial obstacle of passing a member function to the GSL was solved by combining the two approaches explained in the last sections. First of all we define all member function we want to integrate to be static. To still be able, to use attributes from the class, we add a pointer to the struct on which our void pointer points to. Below a minimal working example.

```cpp
#include <iostream>

class my_class {
private:
  double y;

public:
  my_class(double y):y(y) {};
  static double my_func(double x, void* p);
  void print_func(double (*func_pt)(double, void*),double x, void* p);
};

struct Options {
  my_class* class_pt;
};

double my_class::my_func(double x, void* p) {
  Options* opts = static_cast<Options*>(p);
  return x+opts->class_pt->y;
};

void my_class::print_func(double (*func_pt)(double, void*), double x, void* p) {
  std::cout << (*func_pt)(x,p) << std::endl;
};

int main() {
  my_class A(2);
  Options opts;
  opts.class_pt = &A;
  A.print_func(my_class::my_func,3,&opts);
};
```
