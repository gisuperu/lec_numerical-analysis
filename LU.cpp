#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <utility>

namespace lu{
    long double EPSILON = 0.0001; //許容誤差範囲
}

class LU{
private:
    //連立方程式(SimultaneousEquations)の要素
    int variable_amount; //変数数=方程式数
    std::vector<std::vector<long double> > coefficient_matrix; //係数行列(方程式数)*(変数数+1)
public:
    LU(int variable_amount, std::vector<std::vector<long double> > coefficient_matrix);//コンストラクター
    std::vector<std::vector<long double> > copyCoefficientMatrix();
    std::vector<long double> runLU();
    void LUdecomposition(std::vector<std::vector<long double> >& L_matrix, std::vector<std::vector<long double> >& U_matrix);
    void showSimultaneousEquations();
    void showSimultaneousEquations(std::vector<std::vector<long double> > coefficient_matrix);
    void printMatrix(std::vector<std::vector<long double> > matrix);
    void printAnswer(std::vector<long double> answer);
};

//コンストラクター
LU::LU(int variable_amount, std::vector<std::vector<long double> > coefficient_matrix){
    this->variable_amount = variable_amount;
    this->coefficient_matrix = coefficient_matrix;
}

std::vector<std::vector<long double> > LU::copyCoefficientMatrix(){
    std::vector<std::vector<long double> > new_coefficient_matrix;
    for(std::vector<long double> equation : this->coefficient_matrix){
        std::vector<long double> new_equation;
        for(long double coefficient : equation){
            new_equation.push_back(coefficient);
        }
        new_coefficient_matrix.push_back(new_equation);
    }
    return new_coefficient_matrix;
}

/* LU分解法
Ax = b を LUx = b に変形してからxを求める

(1) Ly = bのyを求めてから
(2) Ux = yとしてxを求める

(1) Ly = bのyを求める
Lは下三角行列であるから
y_{i} = (b_{i} - ∑_{k=0~i-1}(l_{i,k} * y_{k}) )/l_{i,i}
で求められる(iは0,1,2, ... ,n-1)

(2) Ux = yとしてxを求める
Uは上三角行列であるから
x_{i} = y_{i} - ∑_{k=i+1~n-1}(u_{i,k} * x_{k})
で求められる(iはn-1,n-2,n-3, ... ,0)
*/
std::vector<long double> LU::runLU(){
    // 与えられた連立方程式を LUx = b とおく.
    std::vector<long double> b_vec(variable_amount, 0);
    std::vector<std::vector<long double> > L_matrix(variable_amount, std::vector<long double>(variable_amount, 0));
    std::vector<std::vector<long double> > U_matrix(variable_amount, std::vector<long double>(variable_amount, 0));
    for(int i = 0; i < variable_amount; i++){
        b_vec.at(i) = coefficient_matrix.at(i).at(variable_amount);
        U_matrix.at(i).at(i) = 1;
    }

    //LU分解
    LU::LUdecomposition(L_matrix, U_matrix);
    printf("L:\n");
    printMatrix(L_matrix);
    printf("U:\n");
    printMatrix(U_matrix);
    
    //＊LU分解の検算
    std::vector<std::vector<long double> > tmp_matrix(variable_amount, std::vector<long double>(variable_amount, 0));
    for(int i = 0; i < variable_amount; i++){
        for(int j = 0; j < variable_amount; j++){
            for(int k = 0; k < variable_amount; k++){
                tmp_matrix.at(i).at(j) += L_matrix.at(i).at(k) * U_matrix.at(k).at(j);
            }
        }
    }
    printf("LU:\n");
    printMatrix(tmp_matrix);

    //L、Uの行列から連立方程式の解を導く
    //(1) Ly = bのyを求める
    std::vector<long double> y_vec(variable_amount, 0);
    for(int i = 0; i < variable_amount; i++){
        //y_{i} = (b_{i} - ∑_{k=0~i-1}(l_{i,k} * y_{k}) )/l_{i,i}
        double long s = 0;
        for(int k = 0; k < i; k++){
            s += L_matrix.at(i).at(k) * y_vec.at(k);
        }
        y_vec.at(i) = (b_vec.at(i) - s) / L_matrix.at(i).at(i);
    }
    //(2) Ux = yとしてxを求める
    std::vector<long double> x_vec(variable_amount, 0);
    for(int i = variable_amount-1; i >= 0; i--){
        //x_{i} = y_{i} - ∑_{k=i+1~n-1}(u_{i,k} * x_{k})
        double long s = 0;
        for(int k = i+1; k < variable_amount; k++){
            s += U_matrix.at(i).at(k) * x_vec.at(k);
        }
        x_vec.at(i) = (y_vec.at(i) - s);
    }

    return x_vec;
}

/*
* 大文字の変数は行列(A,A_{1,1},L,Uなど)
* 小文字の変数はベクトル(a_{1~n-1,0})
* o、Oは全ての要素が0である
A=LUに分解する
A=LUを以下のように解釈できる
[a_{0,0}    , a_{0    ,1~n-1}]   [l_{0,0}    , o_{0    ,1~n-1}][1          , u_{0    ,1~n-1}]
[a_{1~n-1,0}, A_{1~n-1,1~n-1}] = [l_{1~n-1,0}, L_{1~n-1,1~n-1}][o_{1~n-1,0}, U_{1~n-1,1~n-1}]
以上から右辺を計算すると
                                 [l_{0,0}    , l_{0,0}*u_{0    ,1~n-1}                                  ]
                               = [l_{1~n-1,0}, l_{1~n-1,0}*u_{0,1~n-1} + L_{1~n-1,1~n-1}*U_{1~n-1,1~n-1}]
となるため、
(1)  l_{0    ,0    } = a_{0,0}
(2)  l_{1~n-1,0    } = a_{1~n-1,0}
(3)  u_{0    ,1~n-1} = a_{0    ,1~n-1}/l_{0,0}
(4)  A_{1~n-1,1~n-1} = l_{1~n-1,0}*u_{0,1~n-1} + L_{1~n-1,1~n-1}*U_{1~n-1,1~n-1}
と得られる
また(4)より
(4') A' = A_{1~n-1,1~n-1} - l_{1~n-1,0}*u_{0,1~n-1}  とおく
=> 
(5)  A' = L'U' ( = L_{1~n-1,1~n-1}*U_{1~n-1,1~n-1})
次元が下がったLU分解の式(5)が得られるため再起的に分解することで(1)(2)(3)よりL、Uを決定できる。
*/
//LU分解
void LU::LUdecomposition(std::vector<std::vector<long double> >& L_matrix, std::vector<std::vector<long double> >& U_matrix){
    std::vector<std::vector<long double> > new_coefficient_matrix = LU::copyCoefficientMatrix();
    
    //pivotを対角要素上でずらすことで擬似的に次元を下げる
    for(int pivot = 0; pivot < variable_amount; pivot++){
        //(1)  l_{0    ,0    } = a_{0,0}
        L_matrix.at(pivot).at(pivot) = new_coefficient_matrix.at(pivot).at(pivot);

        //(2)  l_{1~n-1,0    } = a_{1~n-1,0}
        for(int i = pivot+1; i < variable_amount; i++){
            L_matrix.at(i).at(pivot) = new_coefficient_matrix.at(i).at(pivot);
        }

        //(3) u_{0    ,1~n-1} = a_{0    ,1~n-1}/l_{0,0}
        for(int i = pivot+1; i < variable_amount; i++){
            U_matrix.at(pivot).at(i) = new_coefficient_matrix.at(pivot).at(i) / L_matrix.at(pivot).at(pivot);
        }

        //(4') A' = A_{1~n-1,1~n-1} - l_{1~n-1,0}*u_{0,1~n-1}
        for(int i = pivot+1; i < variable_amount; i++){
            for(int j = pivot+1; j < variable_amount; j++){
                //l_{1~n-1,0}*u_{0,1~n-1}の(i,j)番地を求める
                long double lu = L_matrix.at(i).at(pivot) * U_matrix.at(pivot).at(j);

                //A - l*u
                new_coefficient_matrix.at(i).at(j) -= lu;
            }
        }
    }
}

void LU::showSimultaneousEquations(){
    printf("連立方程式:\n");
    for(std::vector<long double> equation : this->coefficient_matrix){
        for(int v = 0; v <= variable_amount; v++){
            long double c = equation.at(v);
            if(v == variable_amount){
                printf("\t = ");
            }else if(c >= 0){
                printf("\t + ");
            }else if(c < 0){
                printf("\t - ");
            }
            
            if(v == variable_amount){
                printf("%.6Lf", c);
            }else{
                printf("%.6Lfx_%02d", fabsl(c), v);
            }
        }
        printf("\n");
    }
}
void LU::showSimultaneousEquations(std::vector<std::vector<long double> > coefficient_matrix){
    printf("連立方程式:\n");
    for(std::vector<long double> equation : coefficient_matrix){
        for(int v = 0; v <= variable_amount; v++){
            long double c = equation.at(v);
            if(v == variable_amount){
                printf("\t = ");
            }else if(c >= 0){
                printf("\t + ");
            }else if(c < 0){
                printf("\t - ");
            }
            
            if(v == variable_amount){
                printf("%.6Lf", c);
            }else{
                printf("%.6Lfx_%02d", fabsl(c), v);
            }
        }
        printf("\n");
    }
}

//行列を表示
void LU::printMatrix(std::vector<std::vector<long double> > matrix){
    for(std::vector<long double> rows : matrix){
        for(long double m : rows){
            printf("%.6Lf ", m);
        }
        printf("\n");
    }
}

void LU::printAnswer(std::vector<long double> answer){
    printf("解:\n");
    for(int i = 0; i < variable_amount; i++){
        printf("\tx_%02d = %.6Lf\n", i, answer.at(i));
    }
}


int main(void){
    //連立方程式の定義
    int variable_amount = 3;
    std::vector<std::vector<long double> > coefficient_matrix = {//連立方程式の拡大係数行列
        { 3,  2,  1,  10},
        { 1,  4,  1,  12},
        { 2,  2,  5,  21}
    };
    
    //関数作成
    LU simultaneous_equations(variable_amount, coefficient_matrix);
    simultaneous_equations.showSimultaneousEquations();
    //ガウスジョルダン法の実行
    std::vector<long double> answer = simultaneous_equations.runLU();
    simultaneous_equations.printAnswer(answer);
    return 0;
}