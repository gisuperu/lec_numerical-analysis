from matplotlib import pyplot as plt
import numpy as np

def least_square_method(regression_division, pair_list):
    terms = (regression_division + 1)
    formula = [0 for _ in range(terms)] # 近似曲線(出力)
    matrix = [[] for _ in range(terms)] # 正規方程式

    # 正規方程式生成
    # 左辺(未知変数，近似曲線の係数がある側)の係数計算
    for i in range((terms)*2):
        # 左辺の係数∑x_{k}^{j}の計算(j = 0 ... 2*(近似曲線の係数数))
        sum = 0
        for point in pair_list:
            sum += pow(point[0], i)
        # 係数を行列の対応位置に挿入
        for j in range(max(i-regression_division, 0), min(i+1, terms)):
            matrix[j].append(sum)

    # 右辺(未知変数がない側)の値計算
    sum_list = []
    for i in range(terms):
        sum = 0
        for point in pair_list:
            sum += point[1]*pow(point[0], i)
        sum_list.append(sum)
    for i in range(terms):
        matrix[i].append(sum_list[i])

    formula = gauss_jordan_method(matrix, terms)
    return formula

def gauss_jordan_method(augmented_matrix, variable_amount):
    for i in range(variable_amount):
        # 係数を1に揃える
        for j in range(i, variable_amount):
            c = augmented_matrix[j][i]
            for k in range(i, variable_amount+1):
                augmented_matrix[j][k] /= c
        # 行同士を引く
        for j in range(i+1, variable_amount):
            for k in range(i, variable_amount+1):
                augmented_matrix[j][k] -= augmented_matrix[i][k]
        
    answer = [i[variable_amount] for i in augmented_matrix]
    for i in reversed(range(variable_amount)):
        for j in reversed(range(i+1, variable_amount)):
            answer[i] -= answer[j]*augmented_matrix[i][j]
    return answer

def calc_polynomial(ascending_order, x):
    ans = 0
    for i in range(len(ascending_order)):
        ans += ascending_order[i]*(x**i)
    return ans

def print_polynomial(ascending_order, x, func="y(x)", variable="x"):
    print("{} =".format(func), end="")
    for i in range(len(ascending_order)):
        print(" {0:+}".format(ascending_order[i]), end="")
        
        if i == 1:
            print("{0}".format(variable), end="")
        elif i > 1:
            print("{0}^{{{1}}}".format(variable, i), end="")
    print()

def main():
    # 入力の点集合
    points = tuple((
        (0.5, 10.01),
        (1.0,  8.71),
        (1.5,  7.41),
        (2.0,  6.92),
        (2.5,  5.94)
    ))

    # 最小二乗法の実行(1次式と2次式)
    liner_function = least_square_method(1, points)
    quadratic_function = least_square_method(2, points)

    print("1次式: ", end="")
    print_polynomial(liner_function, 1)
    print("2次式: ", end="")
    print_polynomial(quadratic_function, 2)

    # グラフの描画
    LEFT_RANGE = 0
    RIGHT_RANGE = 5
    # グラフ用に関数を加工
    p = np.linspace(LEFT_RANGE, RIGHT_RANGE, 1000)
    liner_function = calc_polynomial(liner_function, p)
    quadratic_function = calc_polynomial(quadratic_function, p)
    # グラフの設定
    plt.grid(True) #補助線の描画
    plt.xlim(LEFT_RANGE, RIGHT_RANGE) #X軸の範囲を制限
    plt.xlabel("x-axis")
    plt.ylabel("y-axis")

    plt.plot(p, liner_function, label="直線の式", color="g", linestyle="solid")
    plt.plot(p, quadratic_function, label="2次の式", color="r", linestyle="dashed")
    plt.scatter([p[0] for p in points], [p[1] for p in points], color="b")

    plt.legend(prop={"family":"Hiragino sans"}) #凡例の追加
    plt.show()
    plt.close()

if __name__ == '__main__':
    main()