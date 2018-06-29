#!/usr/bin/env python3
# coding: UTF-8

import itertools
import random

#----------------------------↓solver_greedy.py------------------------------

import sys
import math

from common import print_tour, read_input


def distance(city1, city2):
    return math.sqrt((city1[0] - city2[0]) ** 2 + (city1[1] - city2[1]) ** 2)


def solve(cities):
    N = len(cities)

    dist = [[0] * N for i in range(N)]
    for i in range(N):
        for j in range(i, N):
            dist[i][j] = dist[j][i] = distance(cities[i], cities[j])

    current_city = 0
    unvisited_cities = set(range(1, N))
    tour = [current_city]

    while unvisited_cities:
        next_city = min(unvisited_cities,
                        key=lambda city: dist[current_city][city])
        unvisited_cities.remove(next_city)
        tour.append(next_city)
        current_city = next_city
    print("-----print greedyTour-----")
    print(tour)
    return tour

#----------------------------↑solver_greedy.py------------------------------

def makeTour(cities): # 道順の初期値をランダムに(テキトーに決めてくれる)
    firstTour=list(range(len(cities)))
    random.shuffle(firstTour)
    return firstTour

def calcuDist(cities,tour): # 道順を与えると、トータル距離を計算してくれる
    allDist = 0
    for i in range(len(tour)-1):
        allDist += distance(cities[tour[i]], cities[tour[i+1]])
    allDist += distance(cities[tour[0]], cities[tour[len(tour)-1]])
    return allDist


def annealingoptimize(cities, distGreedy, T=100000, cool=0.9999): # hill climb(?) or yakinamasi部分
    # 初期値
    citiesNumber = len(cities)
    citiesNumberIndex = (list(range(0, citiesNumber)))

    count = 1
    while count <= 5:
        tour = makeTour(cities)
        totalDist = calcuDist(cities, tour)
        calculatedTour = tour[:]
        T = 1000
        while T > 0.0001:
            # 値を交換する二つのindexの組み合わせの決め方をinoYakiとは変えてみた
            # やっていることは、ランダムに一点を選んで、その一点のある程度そばにある点の中からもう一点選んで交換してみるという感じ
            # 個人的にはこっちの方が焼きなまし法っぽくっていいのかなあと思ったんだけど実際どうなんだろう
            # index: 選ばれた道順内での周り順の通し番号
            choicedCombi = random.sample(citiesNumberIndex, 1)
            index0 = choicedCombi[0]
            if index0 >= citiesNumber//4 and index0 < citiesNumber-citiesNumber//4:
                indexs = list(range(index0-citiesNumber//4, index0+citiesNumber//4))
            elif index0 < citiesNumber//4:
                indexs = list(range(0,index0+citiesNumber//4))+list(range(index0-citiesNumber//4+citiesNumber+1,citiesNumber))
            else:
                indexs = list(range(index0-citiesNumber//4+1,citiesNumber))+list(range(0,index0+citiesNumber//4-citiesNumber))
            index1 = random.sample(indexs, 1)
            
            # 同じ2点を選んでしまった場合はやり直す
            if index0 != index1[0]:
                # 選ばれた2点を交換
                calculatedTour[index0], calculatedTour[index1[0]] = calculatedTour[index1[0]], calculatedTour[index0]

                # このcalculatedTourがテキトーに二点のcityを入れ替えた後の道順
                newTotalDist = calcuDist(cities, calculatedTour)
                # ↓これの#消すと焼きなましに(?)、pの決め方テキトーです、ググってテキトーに決めた
                p = pow(math.e, -abs(newTotalDist-totalDist)/T) # 温度から確率を定義する

                if newTotalDist < totalDist or random.random() < p: # newTotalDistが小さければ採用する、大きい場合は確率的に採用する
                    print(count, "回目:", "Annealing", totalDist)
                    tour = calculatedTour
                    totalDist = newTotalDist
                else:
                    calculatedTour[index0], calculatedTour[index1[0]] = calculatedTour[index1[0]], calculatedTour[index0]
                T = T * cool # 温度を下げ

        # 3-opt
        if citiesNumber-5 > 0:
            forSaiki = 0
            while forSaiki < citiesNumber*3000:
                trio = []
                citiesNumberIndex = (list(range(0, citiesNumber-5)))
                choiced = random.sample(citiesNumberIndex, 1)
                i = choiced[0]
                citiesNumberIndex = (list(range(i+2, citiesNumber-3)))
                choiced = random.sample(citiesNumberIndex, 1)
                j = choiced[0]
                citiesNumberIndex = (list(range(j+2, citiesNumber-1)))
                choiced = random.sample(citiesNumberIndex, 1)
                k = choiced[0]

                after = [[] for i in range(4)]
                before = distance(cities[tour[i]], cities[tour[i+1]]) + distance(cities[tour[j]], cities[tour[j+1]]) + distance(cities[tour[k]], cities[tour[k+1]])

                # k+1...i-k...j+1-i+1...j 0
                # k+1...i-j+1...k-j...i+1 1
                # k+1...i-j+1...k-i+1...j 2
                # k+1...i-j...i+1-k...j+1 3
                after[0] = [distance(cities[tour[i]], cities[tour[k]]) + distance(cities[tour[j+1]], cities[tour[i+1]]) + distance(cities[tour[k+1]], cities[tour[j]]), 0]
                after[1] = [distance(cities[tour[j+1]], cities[tour[i]]) + distance(cities[tour[k]], cities[tour[j]]) + distance(cities[tour[k+1]], cities[tour[i+1]]), 1]
                after[2] = [distance(cities[tour[j]], cities[tour[k+1]]) + distance(cities[tour[i]], cities[tour[j+1]]) + distance(cities[tour[k]], cities[tour[i+1]]), 2]
                after[3] = [distance(cities[tour[i]], cities[tour[j]]) + distance(cities[tour[j+1]], cities[tour[k+1]]) + distance(cities[tour[k]], cities[tour[i+1]]), 3]
                after.sort()
                min_after, min_after_num = after[0]

                if before > min_after:
                    if min_after_num == 0:
                        connection = [k, i+1]
                    elif min_after_num == 1:
                        connection = [j+1, j]
                    elif min_after_num == 2:
                        connection = [j+1, i+1]
                    elif min_after_num == 3:
                        connection = [j, k]

                    calculatedTour = tour[:i+1]
                    if connection[0] == j+1:
                        calculatedTour.extend(tour[j+1:k+1])
                    elif connection[0] == j:
                        calculatedTour.extend(reversed(tour[i+1:j+1]))
                    elif connection[0] == k:
                        calculatedTour.extend(reversed(tour[j+1:k+1]))

                    if connection[1] == i+1:
                        calculatedTour.extend(tour[i+1:j+1])
                    elif connection[1] == j:
                        calculatedTour.extend(reversed(tour[i+1:j+1]))
                    elif connection[1] == k:
                        calculatedTour.extend(reversed(tour[j+1:k+1]))

                    calculatedTour.extend(tour[k+1:])

                    # print(calculatedTour)
                    newTotalDist = calcuDist(cities, calculatedTour)
                    tour = calculatedTour
                    totalDist = newTotalDist
                    print(count, "回目:", "3 opt", totalDist)
                forSaiki += 1

             # 2-opt
        forSaiki = 0
        while forSaiki < citiesNumber*3000:
            citiesNumberIndex = (list(range(0, citiesNumber-3)))
            choicedCombi = random.sample(citiesNumberIndex, 1)
            index0 = choicedCombi[0]
            citiesNumberIndex = (list(range(index0+2, citiesNumber-1)))
            choicedCombi1 = random.sample(citiesNumberIndex, 1)
            index1 = choicedCombi1[0]

            before = distance(cities[tour[index0]], cities[tour[index0+1]]) + distance(cities[tour[index1]], cities[tour[index1+1]])
            after = distance(cities[tour[index0]], cities[tour[index1]]) + distance(cities[tour[index0+1]], cities[tour[index1+1]])
            if before > after:
                calculatedTour = tour[:index0+1]
                calculatedTour.extend(reversed(tour[index0+1:index1+1]))
                calculatedTour.extend(tour[index1+1:])
                # print(calculatedTour)
                newTotalDist = calcuDist(cities, calculatedTour)
                tour = calculatedTour
                totalDist = newTotalDist
                print(count, "回目:", "2 opt", totalDist)
            forSaiki += 1

        if count == 1:
            bestTour = tour[:]
            bestDist = totalDist
        else:
            if bestDist>totalDist:
                bestTour = tour[:]
                bestDist = totalDist
        print(count, "回目:", totalDist)

        count += 1
    tour = bestTour
    totalDist = bestDist


    if totalDist < distGreedy: # Greedyより結果が良かったら終了する
        print("better than greedy!")
        print("--------the best tour by hill climb---------")
        print(tour)
        print("-------print totalDist--------")
        print(totalDist)
    else:
        print("worse than greedy...")
        print("--------the best tour by hill climb---------")
        print(tour)
        print("-------print totalDist--------")
        print(totalDist)


#----------------------------↓forMain ------------------------------

if __name__ == '__main__':
    # assert len(sys.argv) > 1
    # tour = solve(read_input(sys.argv[1]))
    # print_tour(tour)
    # random.seed(0)
    assert len(sys.argv) > 1
    cities = read_input(sys.argv[1])
    tourGreedy = solve(cities)
    distGreedy = calcuDist(cities, tourGreedy)
    print("-----print distGreedy------")
    print(distGreedy) # ここまでgreedyの実行(別にgreedyの実行は必要ないです、greedyによる算出結果が欲しかっただけ)

    annealingoptimize(cities, distGreedy)

#----------------------------↑forMain------------------------------


#-----------------↓コメント-------------------

#inoYaki2からの変更点は焼きなましの一点選んだ後のランダムな一点を選ぶ範囲の所
#前にcitiesNumber//4 、後ろにcitiesNumber//4 の範囲から選ぶようにした！
#その最初の一点が経路の中でどのような場所にあったとしても

#-----------------↑コメント-------------------


