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
        print(forSaiki)
        T=100000
        while T>0.0001:
            #値を交換する二つのindexの組み合わせの決め方をinoYakiとは変えてみた
            #やっていることは、　ランダムに一点を選んで、その一点のある程度そばにある点の中からもう一点選んで交換してみるという感じ
            #個人的にはこっちの方が焼きなまし法っぽくっていいのかなあと思ったんだけど実際どうなんだろう
            citiesNumber=len(cities)
            citiesNumberIndex=(list(range(0,citiesNumber)))
            choicedCombi=random.sample(citiesNumberIndex,1)
            index0=choicedCombi[0]
            if index0>=citiesNumber//4 and index0<citiesNumber-citiesNumber//4:
                indexs=range(index0-citiesNumber//4,index0+citiesNumber//4)
            elif index0<citiesNumber//4:
                indexs=range(index0+1,index0+1+citiesNumber//2)
            else:
                indexs=range(index0-citiesNumber//2,index0)
            index1=random.sample(indexs,1)
            a=tour[index0] #選ばれたindexのcity
            b=tour[index1[0]] #選ばれたindexのcity2
            calculatedTour=[]

            for city in tour:
                calculatedTour.append(city)
            
            calculatedTour[index1[0]]=a
            calculatedTour[index0]=b
            #このcalculatedTourがテキトーに二点のcityを入れ替えた後の道順
            newTotalDist=calcuDist(cities,calculatedTour)
            #↓これの#消すと焼きなましに(?)、pの決め方テキトーです、ググってテキトーに決めた
            p= pow(math.e, -abs(newTotalDist-totalDist)/T)

            if newTotalDist<totalDist or random.random()<p: #←これの#消すと焼きなましに(?)
                print(totalDist)
                tour=calculatedTour
                totalDist=newTotalDist

            T=T*cool
        forSaiki=0
        while forSaiki<10000:
            citiesNumber=len(cities)
            citiesNumberIndex=(list(range(0,citiesNumber-3)))
            choicedCombi=random.sample(citiesNumberIndex,1)
            index0=choicedCombi[0]
            citiesNumberIndex=(list(range(index0+2,citiesNumber-1)))
            choicedCombi1=random.sample(citiesNumberIndex,1)
            index1=choicedCombi1[0]
            
            before=distance(cities[tour[index0]],cities[tour[index0+1]])+distance(cities[tour[index1]],cities[tour[index1+1]])
            after=distance(cities[tour[index0]],cities[tour[index1]])+distance(cities[tour[index0+1]],cities[tour[index1+1]])
            if before>after:
                calculatedTour=[]
                for city in tour[:index0+1]:
                    calculatedTour.append(city)
                for city in reversed(tour[index0+1:index1+1]):
                    calculatedTour.append(city)
                #print(calculatedTour)
                for city in tour[index1+1:]:
                    calculatedTour.append(city)
                newTotalDist=calcuDist(cities,calculatedTour)
                tour=calculatedTour
                totalDist=newTotalDist
                print(totalDist)
            forSaiki+=1
            
        
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
    random.seed(0)
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


#-----------------↓コメント-------------------
#よくわからない....
#なんか16以上がかなり重くなる...
#何回か実行すると割と良さげな値が出る(?)
#多分、greedyよりよくなかったらやり直すってとこが重そう
#月以降、あんまり時間取れないからあげちゃった
#勝手にやってごめん...
#何なら無視しちゃって大丈夫だよー！！！

#-----------------↑コメント-------------------

#------------------↓問題点----------------------
#初期値も毎回の山登りのためにランダムに更新しているけど、その意味あるのかわからない
#Tの値はもっと大きい方がいいのか、小さい方がいいのかわからない
#コードに無駄があると思われる
#＊hill climbの現状
#何回か "$ python inoYaki.py input_n.csv" を実行しないと、最短?(sample/saのoutput)が得られない(n=0,1,2)
#何回か $ python inoYaki.py input_n.csv" を実行しても、最短?(sample/saのoutput)が得られない(n=3)
# "$ python inoYaki.py input_n.csv" を実行すると謎のエラーが出る(n>=4)

#------------------↑問題点----------------------
