#include "path10.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <iomanip>

/*
* 1. При ошибке завершать выполнение программы
* 2. Потоки. Не используем по возможности. Толкьо print и scan
* 3. Проверка нуля с малым числом eps=e-?
* 4. Избегать проверок на ноль???
* 5. Фулл стеки, не вектора
* 6. DI > 0. Если нет, пользователь дурак.
* 7. Матрица симметрична+
* ТЕСТ:
* Погрешность, k
* LU разложение
* Невязка?
* 1. погр y 100% x полу
* y-x
* ||y-x||/||y|| относ
* b-Ax
* ||b-Ax||/||b|| неяв?
*/
namespace F{
   typedef float real;
   typedef float realscal;
   constexpr const char *fmt = "%10.6E\n";
}
namespace D{
   typedef double real;
   typedef double realscal;
   constexpr const char *fmt = "%18.14E\n";
}
namespace FD{
   typedef float real;
   typedef double realscal;
   constexpr const char *fmt = "%10.6E\n";
}

using namespace std;
using namespace D;
 
int dim = 0;
int profiles = 0;
int n = 0;
int space = 18;
int precision = 14;

void inputDim(string path){
   ifstream inFile;
   inFile.open(path);
   if(inFile){
      inFile>>dim;
      inFile.close();
   } else{
      cout<<"Can't open file:"+path<<endl;
      exit(EXIT_FAILURE);
   }
}
void inputIA(string path, int *arr){
   ifstream inFile;
   inFile.open(path);
   if(inFile){
      for(int i = 0; i<dim+1; i++){
         inFile>>arr[i];
         profiles = arr[i];
      }
      profiles--;
      inFile.close();
   } else{
      cout<<"Can't open file:"+path<<endl;
      //Добавить выход из программы
   }
}
void inputData(string path, real *arr, int size){
   ifstream inFile;
   inFile.open(path);
   if(inFile){
      for(int i = 0; i<size; i++){
         inFile>>arr[i];
      }
      inFile.close();
   } else{
      cout<<"Can't open file:"+path<<endl;
      //Добавить выход из программы
   }
}
void outputX(real *arr){
   FILE *file;
   if(n==0)fopen_s(&file, "X.txt", "w");
   else fopen_s(&file, "X.txt", "a");
   if(file){
   fprintf(file, "---   n = %d   ---\n", n);
   for(int i = 0; i<dim; i++){
      fprintf(file, fmt, arr[i]);
   }
   fprintf(file, "\n\n");
   fclose(file);
   }
   n++;
}
int getPos(int i, int j, int *IA){
   if(j > i) swap(i, j);
   int count = IA[i+1]-IA[i];
   int firstCol = i - count;
   return IA[i]-1 + (j - firstCol);
}
real getElem(int i, int j, int *IA, real *AL, real *AU){
   int count = 0;
   int firstPos = 0;
   real buf = 0;
   if(j > i){
      count = IA[j+1]-IA[j];
      firstPos = j - count;
      if(i<firstPos)return 0;
      buf = AU[IA[j]-1+(i-firstPos)];
   } else{
      count = IA[i+1]-IA[i];
      firstPos = i - count;
      if(j < firstPos)return 0;
      buf = AL[IA[i]-1+(j-firstPos)];
   }
   return buf;
}
void printData(real *DI, real *AL, real *AU, int *IA, real *vec){
   cout<<"\n\nIA:\n";
   for(int i = 0; i<dim+1; i++){
      cout<<IA[i]<<"\n";
   }
   cout<<"\nAL:\n";
   for(int i = 0; i<profiles; i++){
      cout<<AL[i]<<"\n";
   }
   cout<<"\nAU:\n";
   for(int i = 0; i<profiles; i++){
      cout<<AU[i]<<"\n";
   }
   cout<<"\nDI:\n";
   for(int i = 0; i<dim; i++){
      cout<<DI[i]<<"\n";
   }
   cout<<"\nvec:\n";
   for(int i = 0; i<dim; i++){
      cout<<vec[i]<<"\n";
   }
}
void LU(real *DI, real *AL, real *AU, int *IA){
   realscal Utmp = .0;
   realscal Ltmp = .0;
   int index = 0;
   int icount = 0;
   int jcount = 0;
   for(int i = 0; i<dim; i++){
      //DI
      icount = IA[i+1]-IA[i];
      //DI sum
      for(int m = 0; m<icount; m++){
         Utmp += AL[IA[i]-1+m]*AU[IA[i]-1+m];
      }
      Utmp = DI[i]-Utmp;
      DI[i] = sqrt(Utmp);
      Utmp = .0;
      //AL AU
      for(int j = i+1; j<dim; j++){
         jcount = IA[j+1]-IA[j];
         //AL AU sum
         if(i>=j-jcount){
            index = IA[j]+(i-(j-jcount))-1;
            for(int m = 0; m<i; m++){
               if(m>=j-jcount&&m>=i-icount){
                  Utmp += AU[IA[j]-1+m-(j-jcount)]*AL[IA[i]-1+m-(i-icount)];
                  Ltmp += AU[IA[i]-1+m-(i-icount)]*AL[IA[j]-1+m-(j-jcount)];
               }
            }
            //Меняем AU[index] AL[index]
            Utmp = AU[index]-Utmp;
            Ltmp = AL[index]-Ltmp;
            Utmp /= DI[i];
            Ltmp /= DI[i];
            AU[index] = Utmp;
            AL[index] = Ltmp;
            Utmp = .0;
            Ltmp = .0;
         }
      }
   }
}
void Y(real *DI, real *AL, real *AU, int *IA, real *vec){
   int count = 0;
   realscal buf = .0;
   for(int i = 0; i<dim; i++){
      count = IA[i+1]-IA[i];
      //sum
      for(int j = 1; j<=count; j++){
         buf += AL[IA[i+1]-1-j]*vec[i-j];
      }
      buf = vec[i]-buf;
      buf /= DI[i];
      vec[i] = buf;
      buf = .0;
   }
}
void X(real *DI, real *AL, real *AU, int *IA, real *vec){
   int count = 0;
   realscal buf = .0;
   for(int i = dim-1; i>=0; i--){
      for(int j = i+1; j<dim; j++){
         count = IA[j+1]-IA[j];
         if(count!=0 && i>=j-count){
            buf += AU[IA[j+1]-1-(j-i)]*vec[j];
         }
      }
      buf = vec[i]-buf;
      buf /= DI[i];
      vec[i] = buf;
      buf = .0;
   }
}

void generateHilbert(){
   ofstream DIf(path::DI);
   ofstream IAf(path::IA);
   ofstream ALf(path::AL);
   ofstream AUf(path::AU);
   realscal sum = 0;
   int index = 0;
   DIf << 1;
   IAf << 1 << endl << 1;
   if(dim>1){
      index = 1;
   }
   for(int i = 1; i<dim; i++){
      DIf << endl << pow(2*i+1, -1);
      for(int j = 0; j<i; j++){
         sum = pow(i+j+1, -1);
         ALf << sum << endl;
         AUf << sum << endl;
         index++;
      }
      IAf << endl << index;
   }
   profiles = index-1;
   IAf.close();
   ALf.close();
   AUf.close();
   DIf.close();
}

void debugPrint(real a[10][11]){
   for(int i = 0; i<10; i++){
      ostringstream line;
      line << fixed << setprecision(precision);
      for(int j = 0; j<11; j++){
         line << setw(space) << a[i][j];
      }
      line << endl;
      printf("%s", line.str().c_str());
   }
   printf("\n\n");
}

void Gauss(real *DI, real *AL, real *AU, int *IA, real *vec){
   real a[10][11]{};
   for(int i = 0; i<dim; i++){
      for(int j = 0; j<dim; j++){
         if(i==j){
            a[i][j] = DI[i];
         } else{
            if(i<j){
               if(IA[j+1]-IA[j] != 0 && j-(IA[j+1]-IA[j]) <= i){
                  a[i][j] = AU[IA[j+1]-1-(j-i)];
               }
            } else{
               if(IA[i+1]-IA[i] != 0 && i-(IA[i+1]-IA[i]) <= j){
                  a[i][j] = AL[IA[i+1]-1-(i-j)];
               }
            }
         }
      }
      a[i][dim] = vec[i];
   }
   debugPrint(a);

   for(int i = 0; i < dim; i++){
      // Поиск строки с максимальным элементом в текущем столбце
      int maxRow = i;
      for(int k = i + 1; k < dim; k++){
         if(fabs(a[k][i]) > fabs(a[maxRow][i])){
            maxRow = k;
         }
      }

      // Меняем строки местами
      if(maxRow != i){
         for(int k = 0; k <= dim; k++){
            swap(a[i][k], a[maxRow][k]);
         }
      }

      // Проверка на нулевой ведущий элемент
      if(fabs(a[i][i]) < 1e-18){
         cout << "Система не имеет единственного решения.\t"<<n<<endl;
      }

      // Нормализация и исключение
      for(int k = i + 1; k < dim; k++){
         double factor = a[k][i] / a[i][i];
         for(int j = i; j <= dim; j++){
            a[k][j] -= factor * a[i][j];
         }
      }
   }

   // Обратный ход
   for(int i = dim - 1; i >= 0; i--){
      for(int j = i + 1; j < dim; j++){
         a[i][dim] -= a[i][j] * a[j][dim];
      }
      a[i][dim] /= a[i][i];
   }
   for(int i = 0; i<dim; i++){
      vec[i] = a[i][dim];
   }
}

void prog1(){
   inputDim(path::Dim);
   int *IA = new int[dim];
   inputIA(path::IA, IA);
   real *AL = new real[profiles];
   real *AU = new real[profiles];
   real *DI = new real[dim];
   real *vec = new real[dim];
   inputData(path::DI, DI, dim);
   inputData(path::AL, AL, profiles);
   inputData(path::AU, AU, profiles);
   inputData(path::F, vec, dim);
   LU(DI, AL, AU, IA);
   Y(DI, AL, AU, IA, vec);
   X(DI, AL, AU, IA, vec);
   outputX(vec);
}

void vecrog2(){
   generateHilbert();
   prog1();
}

void prog3(){
   inputDim(path::Dim);
   int *IA = new int[dim];
   inputIA(path::IA, IA);
   real *AL = new real[profiles];
   real *AU = new real[profiles];
   real *DI = new real[dim];
   real *vec = new real[dim];
   inputData(path::DI, DI, dim);
   inputData(path::AL, AL, profiles);
   inputData(path::AU, AU, profiles);
   inputData(path::X, vec, dim);
   printData(DI, AL, AU, IA, vec);
}

void PrintDenseMatrix(real *DI, real *AL, real *AU, int *IA){
   FILE *file;
   if(fopen_s(&file, "A.txt", "w") != 0){
      cerr << "Ошибка: не удалось открыть файл " << endl;
      return;
   }

   for(int i = 0; i < dim; ++i){
      ostringstream line;
      line << fixed << setprecision(precision);
      for(int j = 0; j < dim; ++j){
         if(i==j){
            line << setw(space) << DI[i];
         } else{
            line << setw(space) << getElem(i,j,IA,AL,AU);
         }
      }
      line << endl;
      fprintf(file, "%s", line.str().c_str());
   }
   fclose(file);
}

void testmaker(){
   inputDim(path::Dim);
   int *IA = new int[dim];
   inputIA(path::IA, IA);
   real *AL = new real[profiles];
   real *AU = new real[profiles];
   real *DI = new real[dim];
   real *vec = new real[dim];

   real *d1 = new real[19];
   real *f1 = new real[19];
   inputData("d1.txt", d1, 19);
   inputData("f1.txt", f1, 19);

   for(int i = 0; i<1; i++){
      inputData(path::DI, DI, dim);
      inputData(path::AL, AL, profiles);
      inputData(path::AL, AU, profiles);
      inputData(path::F, vec, dim);
      DI[0] = d1[i];
      vec[0] = f1[i];
      //LU(DI, AL, AU, IA);
      //Y(DI,AL,AU,IA,vec);
      //X(DI,AL,AU,IA,vec);
      Gauss(DI,AL,AU,IA,vec);
      outputX(vec);
   }
}

int main(){
   testmaker();
   return 0;
}