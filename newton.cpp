#include<stdio.h>
#include<iostream>
#include<vector>

namespace newton{
    const int ROOT = 0; // 関数の解(X軸の位置)
    const long double EPSILON = 0.0001; //許容誤差範囲
    const int ETERNAL_ROOP_LIMIT = 100000000; //再帰還数(ニュートン法の実行メソッド)の最大ループ回数
}


class Newton{
private:
    int division; //関数の次元
    std::vector<int> coefficients; //係数行列(昇べきの順)
public:
    void set(int division, std::vector<int> coefficients);
    void printFunction(); //関数表示用メソッド
    long double function(long double x);
    long double differentiatedFunction(long double x);
    long double run(long double a, int roop_counter);
};

//コンストラクター
void Newton::set(int division, std::vector<int> coefficients){
    this->division = division;
    this->coefficients = coefficients;
}

//保持している情報から関数の多項式表示
void Newton::printFunction(){
    printf("f(x) = ");
    for(int d = division; d >= 0; d--){
        int c = coefficients.at(d);
        if(c == 0){
            continue;
        }else if(c > 0 && d != division){
            printf(" + ");
        }else if(c < 0 && d != division){
            printf(" - ");
        }
        
        if(c == 1){
            printf("x^%d", d);
        }else if(d == 0){
            printf("%d", abs(c));
        }else{
            printf("%dx^%d", abs(c), d);
        }
    }
    printf("\n");
}

//関数から値を返す
long double Newton::function(long double x){
    long double result = 0;
    long double weight = 1.0;
    for(int c : coefficients){
        result += c * weight;
        weight *= x;
    }
    return result;
}

//導関数の値を返す
long double Newton::differentiatedFunction(long double x){
    long double result = 0;
    long double weight = 1.0;
    int dc = 1; //微分で増える各項の係数値
    for(int i = 1; i <= division; i++){
        dc = i;
        result += coefficients.at(i) * dc * weight;
        weight *= x;
    }
    return result;
}

//Newton法の実行
long double Newton::run(long double a, int roop_counter){
    if(roop_counter > newton::ETERNAL_ROOP_LIMIT){ // 再帰上限
        std::cerr << "error : 再帰回数が上限に達しました" << std::endl;
        return 0;
    }

    long double b = a - Newton::function(a)/Newton::differentiatedFunction(a); //aの接線とx軸との交点
    long double r = a - b; //区間の差
    r = r > 0 ? r : -r; //絶対値
    // debug log(begin) ----
    std::cerr << "(a, b) = " << "(" << a << ", " << b << ")" << std::endl;
    // ---- debug log(end)
    if(r < newton::EPSILON){ //誤差EPSILON以下は終了
        return b;
    }else{
        return Newton::run(b, roop_counter+1);
    }
}

int main(){
    //関数の指定
    int division = 2; //関数の次元
    int org_data[] = {-2, 0, 1}; //初期値用係数行列(昇べきの順)
    std::vector<int> coefficients; //係数行列(昇べきの順)
    for(int i : org_data){
        coefficients.push_back(i);
    }

    // 始点の指定
    long double a = 1000; //区間の下限
    
    //関数作成
    Newton fx; 
    fx.set(division, coefficients);

    fx.printFunction();
    printf("近似解 = %Lf\n", fx.run(a, 0));

    return 0;
}
