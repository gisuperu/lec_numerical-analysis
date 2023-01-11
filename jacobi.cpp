#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <utility>

namespace jacobi{
    long double EPSILON = 0.0001; //許容誤差範囲
    int MAX_LOOP = 50; //最大繰り返し回数
}

class Jacobi{
private:
    //連立方程式(SimultaneousEquations)の要素
    int variable_amount; //変数数=方程式数
    std::vector<std::vector<long double> > coefficient_matrix; //係数行列(方程式数)*(変数数+1)
public:
    Jacobi(int variable_amount, std::vector<std::vector<long double> > coefficient_matrix);//コンストラクター
    std::vector<std::vector<long double> > copyCoefficientMatrix();
    std::vector<long double> runJacobi();
    long double runAjustEquation(std::vector<long double> equation, std::vector<long double> equation_parameter, int variable_number);
    void showSimultaneousEquations();
    void showSimultaneousEquations(std::vector<std::vector<long double> > coefficient_matrix);
    void printAnswer(std::vector<long double> answer);
};

//コンストラクター
Jacobi::Jacobi(int variable_amount, std::vector<std::vector<long double> > coefficient_matrix){
    this->variable_amount = variable_amount;
    this->coefficient_matrix = coefficient_matrix;
}

std::vector<std::vector<long double> > Jacobi::copyCoefficientMatrix(){
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

std::vector<long double> Jacobi::runJacobi(){
    std::vector<std::vector<long double> > new_coefficient_matrix = Jacobi::copyCoefficientMatrix();
    std::vector<long double>               answer                 = {1, 1, 1}; //解の初期値 

    // 修正式を用いて解の計算
    for(int loop = 0; loop < jacobi::MAX_LOOP; loop++){
        std::vector<long double> next_answer(answer);
        for(int i = 0; i < variable_amount; i++){
            std::vector<long double> equation = new_coefficient_matrix.at(i);

            next_answer.at(i) = runAjustEquation(equation, answer, i);//修正式
        }
        std::cout << loop+1 << "回目" << std::endl;
        Jacobi::printAnswer(next_answer);

        // 絶対値誤差の総和
        long double difference = 0;
        for(int i = 0; i < variable_amount; i++){
            difference += fabsl(next_answer.at(i) - answer.at(i));
        }
        
        // 許容誤差範囲なら終了
        if(difference < jacobi::EPSILON){
            return next_answer;
        }

        answer = next_answer;
    }
    printf("最大繰り返し回数を超過しました\n");
    //解の出力
    return answer;
}

//修正式の計算
long double Jacobi::runAjustEquation(std::vector<long double> equation, std::vector<long double> equation_parameter, int variable_number){
    long double answer = equation.at(variable_amount);
    for(int i = 0; i < variable_amount; i++){
        if(i == variable_number){
            continue;
        }
        answer -= equation.at(i)*equation_parameter.at(i);
    }
    answer /= equation.at(variable_number);
    return answer;
}

void Jacobi::showSimultaneousEquations(){
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
void Jacobi::showSimultaneousEquations(std::vector<std::vector<long double> > coefficient_matrix){
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

void Jacobi::printAnswer(std::vector<long double> answer){
    printf("解: ");
    for(int i = 0; i < variable_amount; i++){
        printf("\tx_%02d = %.6Lf", i, answer.at(i));
    }
    printf("\n");
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
    Jacobi simultaneous_equations(variable_amount, coefficient_matrix);
    simultaneous_equations.showSimultaneousEquations();
    //ガウスザイデル法の実行
    std::vector<long double> answer = simultaneous_equations.runJacobi();
    simultaneous_equations.printAnswer(answer);
    return 0;
}