#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <utility>

namespace gaussJordan{
    long double EPSILON = 0.0001; //許容誤差範囲
}

class GaussJordan{
private:
    //連立方程式(SimultaneousEquations)の要素
    int variable_amount; //変数数=方程式数
    std::vector<std::vector<long double> > coefficient_matrix; //係数行列(方程式数)*(変数数+1)
public:
    GaussJordan(int variable_amount, std::vector<std::vector<long double> > coefficient_matrix);//コンストラクター
    std::vector<std::vector<long double> > copyCoefficientMatrix();
    std::vector<long double> runGaussJordan(); 
    void showSimultaneousEquations();
    void showSimultaneousEquations(std::vector<std::vector<long double> > coefficient_matrix);
    void printAnswer(std::vector<long double> answer);
};

//コンストラクター
GaussJordan::GaussJordan(int variable_amount, std::vector<std::vector<long double> > coefficient_matrix){
    this->variable_amount = variable_amount;
    this->coefficient_matrix = coefficient_matrix;
}

std::vector<std::vector<long double> > GaussJordan::copyCoefficientMatrix(){
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

std::vector<long double> GaussJordan::runGaussJordan(){
    std::vector<std::vector<long double> > new_coefficient_matrix = GaussJordan::copyCoefficientMatrix();
    //係数行列の対角要素を1にする
    for(int i = 0; i < variable_amount; i++){
        int pivot = i; 
        int p = i; //pivot番目の項が存在する方程式の行
        for(p = i; p < variable_amount; p++){
            if(new_coefficient_matrix.at(p).at(pivot) != 0){
                break;
            }
        }
        if(p >= variable_amount){//解が一意に決まらない場合
                std::cerr << "解が一意に定まりません" << std::endl;
                std::vector<long double> v;
                return v;
        }
        // std::cerr << "(pivot, p) = (" << pivot << ", " << p << ")" << std::endl;
        std::swap(new_coefficient_matrix.at(pivot), new_coefficient_matrix.at(p));

        //pivot番目の項の係数を1にする．
        for(int j = pivot; j < variable_amount; j++){
            long double pivot_coefficient = new_coefficient_matrix.at(j).at(pivot);
            //pivot番目の項の係数が既に0で存在しない場合はその方程式を飛ばす
            if(pivot_coefficient == 0){
                continue;
            }

            for(int k = pivot; k < variable_amount + 1/*一つの方程式の項数*/; k++){
                new_coefficient_matrix.at(j).at(k) /= pivot_coefficient;
            }
        }

        //pivot行の方程式残して他のpivot番目の係数を0にするように
        //pivot行の方程式と差をとる．
        for(int j = pivot+1/**/; j < variable_amount; j++){
            long double c = new_coefficient_matrix.at(j).at(pivot);
            for(int k = pivot; k < variable_amount + 1/*一つの方程式の項数*/; k++){
                if(c != 0){
                    new_coefficient_matrix.at(j).at(k) -= new_coefficient_matrix.at(pivot).at(k);
                }
            }
        }
        showSimultaneousEquations(new_coefficient_matrix);
    }
    //解のベクトルを出力
    std::vector<long double> answer(variable_amount);
    for(int e = variable_amount-1; e >= 0; e--){
        long double ans = new_coefficient_matrix.at(e).at(variable_amount);
        for(int c = e+1; c < variable_amount; c++){
            ans -= new_coefficient_matrix.at(e).at(c)*answer.at(c);
        }
        answer.at(e) = ans;
    }
    return answer;
}

void GaussJordan::showSimultaneousEquations(){
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
void GaussJordan::showSimultaneousEquations(std::vector<std::vector<long double> > coefficient_matrix){
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

void GaussJordan::printAnswer(std::vector<long double> answer){
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
    GaussJordan simultaneous_equations(variable_amount, coefficient_matrix);
    simultaneous_equations.showSimultaneousEquations();
    //ガウスジョルダン法の実行
    std::vector<long double> answer = simultaneous_equations.runGaussJordan();
    simultaneous_equations.printAnswer(answer);
    return 0;
}