#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
void output_test(NumericMatrix test){
  std::ofstream output;
  string filename = "greb.csv";
  output.open(filename.c_str());

  for(int i = 0; i <test.nrow();++i){
    for(int j = 0; j < test.ncol();++j){
      output << test(i,j) << ",";
    }
    output << endl;
  }
  output.close();
}

//[[Rcpp::export]]
void mod_test(NumericVector test){
  for(int i = 0;i < test.size();++i){
    if(i%2==0){
      cout << test(i) << endl;
    }
  }
}

//[[Rcpp::export]]
NumericMatrix copy_test(NumericMatrix test){
  NumericMatrix final;
  final = test;
  return(final);
}

//[[Rcpp::export]]
void another_test(NumericVector startvalue, double sampno, double probab){
  ofstream csv_write;
  csv_write.open("test.csv");
  csv_write << sampno << ",";
  for(NumericVector::iterator it = startvalue.begin();it != startvalue.end();++it){
    csv_write << *it << ",";
  }
  csv_write << probab << endl;
  csv_write.close();
}
//[[Rcpp::export]]
double rnorm_test(double current, double step){
  return(R::rnorm(current,step));
}


//[[Rcpp::export]]
double get_matrix_index(int row, int column){
  double my_array[3][3];
  my_array[0][0] = 1;
  my_array[0][1] = 2;
  my_array[0][2] = 3;
  my_array[1][0] = 4;
  my_array[1][1] = 5;
  my_array[1][2] = 6;
  my_array[2][0] = 7;
  my_array[2][1] = 8;
  my_array[2][2] = 9;
  cout << "Size: " << sizeof(my_array)/sizeof(my_array[0][0]) << endl;
  return(my_array[row][column]);

}

//[[Rcpp::export]]
int numeric_vector(NumericVector test, int iterations){
  double tmp;
  for(int j =0;j<iterations;++j){
    for(int  i = 0; i < test.size(); ++i){
      tmp = test(i);
    }
  }
  return(test.size());
}

//[[Rcpp::export]]
int numeric_vector2(NumericVector test, int iterations){
  double tmp;
  double test2[test.size()];
  for(int i = 0; i < test.size();++i){
    test2[i] = test(i);
  }
  for(int j =0;j<iterations;++j){
    for(int  i = 0; i < test.size(); ++i){
      tmp = test2[i];
    }
  }
  return(test.size());
}


//[[Rcpp::export]]
int numeric_vector3(NumericVector test, int iterations){
  double tmp;
  for(int j =0;j<iterations;++j){
    for(NumericVector::iterator it = test.begin();it != test.end();++it){
      tmp = *it;
    }
  }
  return(test.size());
}
//[[Rcpp::export]]
int numeric_matrix(NumericMatrix test, int iterations){
  double tmp;
  for(int x = 0; x < iterations;++x){
    for(int i = 0;i<test.nrow();++i){
      for(int j =0; j < test.ncol();++j){
	tmp = test(i,j);
      }
    }
  }
  return(test.nrow());
}

//[[Rcpp::export]]
int numeric_matrix2(NumericMatrix test, int iterations){
  double tmp;
  double test2[test.nrow()][test.ncol()];
  for(int x = 0; x < iterations;++x){
  for(int i = 0;i<test.nrow();++i){
    for(int j =0; j < test.ncol();++j){
      tmp = test2[i][j];
    }
  }
  }
  return(test.nrow());
}


//[[Rcpp::export]]
double my_dnorm(double x, double mean, double sd, bool log){
  return(R::dnorm(x,mean,sd,log));
}
