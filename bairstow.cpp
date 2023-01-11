#include<stdio.h>
#include<iostream>
#include<vector>
#include<cmath>

/* --- ---　Bairstow's methodの概要 --- ---
与えられる関数
f(x) =  a_0 + a_1*x + a_2*x^2 + ... + a_n*x^n
     = (b_0 + b_1*x + b_2*x^2 + ... + b_{n-2}*x^{n-2})*(x^2 + px + q) + Ax + B
とおく。
f(x) = (b_0q + B) + (b_0p + b_1q + A)x + (b_0 + b_1p + b_2q)x^2 + ... + (b_{n-2} + b_{n-1}p + b_nq)x^n
但し、b_{n-1} = b_n = 0  とする。
よって、
=> [1] a_0 = b_0q + B
   [2] a_1 = b_0p + b_1q + A
   [3] a_k = b_{k-2} + b_{k-1}p + b_kq (kはn以下の自然数とし、b_{-1} = b_{-2} = 0とする。(要検証、多分いらない))
以上三式から
=> [4] A = b_{-1}
   [5] B = b_{-2} + b_{-1}p
AとBはpとqの関数であるからテイラー展開より
   [6] A(p_0 + ∆p, q_0 + ∆q) = A(p_0, q_0) + ∂A/∂p *∆p + ∂A/∂q *∆q = 0
   [7] B(p_0 + ∆p, q_0 + ∆q) = B(p_0, q_0) + ∂B/∂p *∆p + ∂B/∂q *∆q = 0
が得られる。(上二式が0になる時が(x^2 + px + q)で因数分解できた場合である。)
[3]より、
   [8]  ∂b_k/∂p = -b_{k+1} - ∂b_{k+1}/∂p *p - ∂b_{k+2}/∂p *q
   [9]  ∂b_k/∂q = -b_{k+2} - ∂b_{k+1}/∂q *p - ∂b_{k+2}/∂q *q
   [10] ∂b_k/∂p = ∂b_{k-1}/∂q  が得られるため
   [11] c_k = -∂b_{k-1}/∂p = ∂b_{k-2}/∂q
とおく。
[8]または[9]、加えて[11]より
   [12] c_k = b_k -p*c_{k+1} -q*c_{k+2}
が得られる。
[4][5]、[10][11][12]より[6][7]を連立方程式として∆p,∆qを解くと、
   [13] ∆p = (c_0*b_{-1} - c_1*b_{-2})               / (c_0^2 + c_1*(b_{-1} - c_{-1}))
   [14] ∆q = (c_0*b_{-2} + b_{-1}*(b_{-1} - c_{-1})) / (c_0^2 + c_1*(b_{-1} - c_{-1}))
が得られる。
(Bairstow's methodでは∆pと∆qが許容誤差(イプシロン)以下になるまで
p_k = p_{k-1} + ∆p, q_k = q_{k-1} + ∆qを繰り返して漸近していく。
求め終われば(b_0 + b_1*x + b_2*x^2 + ... + b_{n-2}*x^{n-2})に対して同様に繰り返して
解の全てを求めていく。)
--- --- --- --- */
namespace bairstow{
    const int ROOT = 0; // 関数の解(X軸の位置)
    const long double EPSILON = 0.0001; //許容誤差範囲
}

class Bairstow
{
private:
    int division; //関数の次元
    std::vector<long double> coefficients; //係数行列(昇べきの順)
    std::vector<std::pair<long double, long double> > answers; //解の配列(pairは実部と虚部)
public:
    void set(int division, std::vector<long double> coefficients);
    void printFunction(); //関数表示用メソッド
    void calcRoot(long double p, long double q); //answersに直接、解が書き込まれる
    void run(); //実行結果はanswersに直接、解が書き込まれる
    std::vector<std::pair<long double, long double> > getAnswers();
};

//コンストラクター
void Bairstow::set(int division, std::vector<long double> coefficients){
    this->division = division;
    this->coefficients = coefficients;
}

//保持している情報から関数の多項式表示
void Bairstow::printFunction(){
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

//二次方程式の解を返す(x^2 + p*x + q)
void Bairstow::calcRoot(long double p, long double q){
    long double D = p*p -4 * q;//判定式D
    long double r1, r2; //解の実部
    long double i1, i2; //解の虚部
    std::cerr << "(p, q) = (" << p << ", " << q << ")" << std::endl;
    std::cerr << "D = " << D << std::endl;
    if(D >= 0){ //実数解(重解は同じ解が二つでる)
        long double sqrt_D = std::sqrt(D);
        r1 = -(p - sqrt_D)/2.0;
        r2 = -(p + sqrt_D)/2.0;
        i1 = 0.0;
        i2 = 0.0;
        this->answers.push_back(std::make_pair(r1, i1));
        this->answers.push_back(std::make_pair(r2, i2));
    }else{
        long double sqrt_D = std::sqrt(-D);
        r1 = -p/2.0;
        r2 = -p/2.0;
        i1 = sqrt_D/2.0;
        i2 = -sqrt_D/2.0;
        this->answers.push_back(std::make_pair(r1, i1));
        this->answers.push_back(std::make_pair(r2, i2));
    }
}

//Bairstow's methodの実行
void Bairstow::run(){
    int next_division = division;
    std::vector<long double> next_coefficients;
    for (long double &c : coefficients) {
        next_coefficients.push_back(c);
    }
    while(true){
        if(next_division <= 0){
            return ;
        }else if(next_division == 1){
            answers.push_back(std::make_pair(-next_coefficients.at(0)/next_coefficients.at(1), 0));
            return ;
        }else if(next_division == 2){
            Bairstow::calcRoot(next_coefficients.at(1)/next_coefficients.at(2), next_coefficients.at(0)/next_coefficients.at(2));
            return ;
        }else{ //next_coefficientsから二次式(x^2 + px + q)を因数として出す
            long double p = 1; //初期値1
            long double q = 1; //初期値1
            int max_item = (next_division + 1) + 2;
            std::vector<long double> b(max_item, 0); //末尾2要素はb_{-1}とb_{-2}を格納するため、インデックスはmax_itemの剰余を用いる
            std::vector<long double> c(max_item, 0); //末尾2要素はc_{-1}とc_{-2}を格納するため、インデックスはmax_itemの剰余を用いる
            
            long double dp, dq;
            do{
                for(int i = next_division-2; i >= -2; i--){ // b_{division-1} = b_{division} = 0を初期値として漸化式でbを埋める
                    b.at((i+max_item)%max_item) = next_coefficients.at((i+2+max_item)%max_item) - p*b.at((i+1+max_item)%max_item) - q*b.at((i+2+max_item)%max_item);
                }
                for(int i = next_division-2; i >= -2; i--){ // c_{division-1} = c_{division} = 0を初期値として漸化式でcを埋める
                    c.at((i+max_item)%max_item) = b.at((i+max_item)%max_item) - p*c.at((i+1+max_item)%max_item) - q*c.at((i+2+max_item)%max_item);
                }
                long double denominator = std::pow(c.at(0), 2) + c.at(1) * (b.at((-1+max_item)%max_item) - c.at((-1+max_item)%max_item));
                dp = (c.at(0) * b.at((-1+max_item)%max_item) - c.at(1) * b.at((-2+max_item)%max_item)) / denominator;
                dq = (c.at(0) * b.at((-2+max_item)%max_item) + b.at((-1+max_item)%max_item) * (b.at((-1+max_item)%max_item) - c.at((-1+max_item)%max_item))) / denominator;
                
                // debug log(begin) ----
                std::cerr << "(dp, dp) = (" << dp << ", " << dq << ")" << std::endl;
                // ---- debug log(end)
                p += dp;
                q += dq;
            }while((std::abs(dp) > bairstow::EPSILON) && (std::abs(dq) > bairstow::EPSILON));
            Bairstow::calcRoot(p, q);
            next_division -= 2;
            std::vector<long double> new_coefficients;
            for(int i = 0; i <= next_division; i++){
                new_coefficients.push_back(b.at(i));
            }
            next_coefficients = new_coefficients;
        }
    }
}

//getter(answers)
std::vector<std::pair<long double, long double> > Bairstow::getAnswers(){
    return this->answers;
}



int main(){
    //関数の指定
    int division = 2; //関数の次元
    long double org_data[] = {-2, 0, 1}; //初期値用係数行列(昇べきの順)
    std::vector<long double> coefficients; //係数行列(昇べきの順)
    for(int i : org_data){
        coefficients.push_back(i);
    }

    
    //関数作成
    Bairstow fx; 
    fx.set(division, coefficients);

    fx.printFunction();
    fx.run();
    std::vector<std::pair<long double, long double> > answers = fx.getAnswers();
    for(std::pair<long double, long double> a : answers){
        if(a.second >= 0){
            printf("近似解 = %Lf + i* %Lf\n", a.first, a.second);
        }else{
            printf("近似解 = %Lf - i* %Lf\n", a.first, -a.second);
        }
    }

    return 0;
}
