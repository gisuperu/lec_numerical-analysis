#include<stdio.h>
#include<iostream>
#include<vector>

namespace nibun{
    const int ROOT = 0; // 関数の解(X軸の位置)
    const long double EPSILON = 0.0001; //許容誤差範囲
}


class Nibun{
private:
    int division; //関数の次元
    std::vector<int> coefficients; //係数行列(昇べきの順)
public:
    void set(int division, std::vector<int> coefficients);
    void printFunction(); //関数表示用メソッド
    long double function(long double x);
    long double run(long double range_lower, long double range_higher);
};

//コンストラクター
void Nibun::set(int division, std::vector<int> coefficients){
    this->division = division;
    this->coefficients = coefficients;
}

//保持している情報から関数の多項式表示
void Nibun::printFunction(){
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
long double Nibun::function(long double x){
    long double result = 0;
    long double weight = 1.0;
    for(int c : coefficients){
        result += c * weight;
        weight *= x;
    }
    return result;
}

//二分法の実行
long double Nibun::run(long double range_lower, long double range_higher){
    long double c = (range_higher + range_lower) / 2; //区間の中点
    long double fa = Nibun::function(range_lower); //区間下限の関数値
    long double fb = Nibun::function(range_higher); //区間上限の関数値
    long double fc = Nibun::function(c); // 中点の関数値
    long double r = range_higher - range_lower; //区間の差
    r = r > 0 ? r : -r; //絶対値
    // debug log(begin) ----
    std::cerr << "(a, c, b) = " << "(" << range_lower << ", " << c << ", " << range_higher << ")" << std::endl;
    // ---- debug log(end)
    if(r < nibun::EPSILON){ //誤差EPSILON以下は終了
        return c;
    }
    if(fc > nibun::ROOT){
        return Nibun::run(range_lower, c);
    }else if(fc < nibun::ROOT){
        return Nibun::run(c, range_higher);
    }else {
        return c;
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

    // 区間の指定
    long double range_lower = 0; //区間の下限
    long double range_higher = 1000; //区間の上限
    
    //関数作成
    Nibun fx; 
    fx.set(division, coefficients);

    fx.printFunction();
    printf("近似解 = %Lf\n", fx.run(range_lower, range_higher));

    return 0;
}
