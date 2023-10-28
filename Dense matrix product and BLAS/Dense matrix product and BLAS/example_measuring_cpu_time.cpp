/// this is to illustrate mesuring cpu time
// to compile you can for example type :
// g++ -std=c++11 example_measuring_cpu_time.cc
// you can then run it typeing ./a.out
#include <chrono>
#include <iostream>

int main(){
  //define two variables to store time points :
  std::chrono::time_point<std::chrono::system_clock> start, end;
  // set start to the current time
  start = std::chrono::system_clock::now();
  // Wait unti the user type return
  char a;
  std::cin >> a;
  // set end to the current time
  end   = std::chrono::system_clock::now();
  // save the time in second between start and end
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "I waited " << elapsed_seconds.count() << " Seconds "<< std::endl;
}
