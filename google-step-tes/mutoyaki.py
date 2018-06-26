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
def makeTour(cities):#道順の初期値をランダムに(テキトーに決めてくれる)
    firstTour=list(range(len(cities)))
    random.shuffle(firstTour)
    return firstTour


def calcuDist(cities,tour):#道順を与えると、トータル距離を計算してくれる
    allDist=0
    for i in range(len(tour)-1):
        allDist+=distance(cities[tour[i]],cities[tour[i+1]])
    allDist+=distance(cities[tour[0]],cities[tour[len(tour)-1]])
    return allDist


def annealingoptimize(cities,firstTour,allDist,distGreedy,T=100000, cool=0.9999):#hill climb(?) or yakinamasi部分
    forSaiki=0
    while forSaiki<100001:
        #初期値
        tour=makeTour(cities)
        totalDist=calcuDist(cities,tour)
        T=10
        while T>0.0001:
            #値を交換する二つのindexの組み合わせをランダムに決める
            #下記テキトーな二つのcityを入れ替えるための操作
            citiesNumber=len(cities)
            citiesNumberIndex=(list(range(0,citiesNumber-3)))
            choicedCombi=random.sample(citiesNumberIndex,1)
            index0=choicedCombi[0]
            citiesNumberIndex=(list(range(index0+2,citiesNumber-1)))
            choicedCombi1=random.sample(citiesNumberIndex,1)
            index1=choicedCombi1[0]
            a=tour[index0] #選ばれたindexのcity
            b=tour[index1] #選ばれたindexのcity2
            print(tour[index0])
            print(tour[index0+1])
            print(tour[index1])
            print(tour[index1+1])
            before=distance(cities[tour[index0]],cities[tour[index0+1]])+distance(cities[tour[index1]],cities[tour[index1+1]])
            after=distance(cities[tour[index0]],cities[tour[index1]])+distance(cities[tour[index0+1]],cities[tour[index1+1]])
            print(tour)
            calculatedTour=[]
            for city in tour[:index0+1]:
                calculatedTour.append(city)
            print(calculatedTour)
            for city in reversed(tour[index0+1:index1+1]):
                calculatedTour.append(city)
            print(calculatedTour)
            for city in tour[index1+1:]:
                    calculatedTour.append(city)
            print(calculatedTour)
            
            #このcalculatedTourがテキトーに二点のcityを入れ替えた後の道順
            newTotalDist=calcuDist(cities,calculatedTour)
            #↓これの#消すと焼きなましに(?)、pの決め方テキトーです、ググってテキトーに決めた
            p= pow(math.e, -abs(newTotalDist-totalDist)/T)

            if after<before or random.random()<p: #←これの#消すと焼きなましに(?)
                print(totalDist)
                tour=calculatedTour
                totalDist=newTotalDist

            T=T*cool
        
        if totalDist<distGreedy:#Greedyより結果が良かったら終了する
            print("--------the best tour by hill climb---------")
            print(tour)
            print("-------print totalDist--------")
            print(totalDist)
            break
        forSaiki+=1
        if forSaiki==100000:
            print("break")
    
        

#----------------------------↓forMain ------------------------------
if __name__ == '__main__':
    #assert len(sys.argv) > 1
    #tour = solve(read_input(sys.argv[1]))
    #print_tour(tour)
    assert len(sys.argv) > 1
    cities=read_input(sys.argv[1])
    tourGreedy = solve(cities)
    distGreedy=calcuDist(cities,tourGreedy)
    print("-----print distGreedy------")
    print(distGreedy)#ここまでgreedyの実行(別にgreedyの実行は必要ないです、greedyによる算出結果が欲しかっただけ)

    firstTour=makeTour(cities)
    firstDist=calcuDist(cities,firstTour)
    annealingoptimize(cities,firstTour,firstDist,distGreedy)
   
#----------------------------↑forMain------------------------------




#------------------↓問題点----------------------


#------------------↑問題点----------------------
