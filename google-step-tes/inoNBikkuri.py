import csv
import itertools
import math

def csvToXyList():
    csv_file = open("input_0.csv","r", encoding="ms932",errors="",newline="")#csvファイルは適宜変える
    f=csv.reader(csv_file, delimiter="\t", doublequote=True, lineterminator="\r\n", quotechar='"',skipinitialspace=True)
    header=next(f)
    for row in f:
        rowlist=row[0].split(",")
        x.append(float(rowlist[0]))
        y.append(float(rowlist[1]))
    # x=[x1,x2,x3,...]
    # y=[y1,y2,y3,...]

def makeIndexs():
    dataLength=len(x)#=len(y)
    indexs=range(0,dataLength)#indexs=[0,1,2,3...dataLength-1]
    return indexs

def SearchMinimumDist(indexs):
    dataLength=len(x)#=len(y)
    minimumCombi=[]
    minimumDis=0
    for combi in itertools.permutations(indexs):#すべての組み合わせを試す
        #combi=[4,2,0,3,1] (ex) 
        #combi[0]=4 (ex)
        #combi[1]=2 (ex)
        #combi[2]=0 (ex)
        #combi[3]=3 (ex)
        #combi[4]=1 (ex)
        #dis1=(x[combi[0]]-x[combi[1]])*(x[combi[0]]-x[combi[1]])+(y[combi[0]]-y[combi[1]])*(y[combi[0]]-y[combi[1]])

        #for distance
        param=0
        disz=[]#それぞれの二点間の距離を入れる用のリスト
        while param < dataLength-1:
            disz.append(math.sqrt((x[combi[param]]-x[combi[param+1]])**2+(y[combi[param]]-y[combi[param+1]])**2))
            param+=1
        #一番最後と一番最初の点の距離だけ特別に計算
        disLast= math.sqrt((x[combi[0]]-x[combi[dataLength-1]])**2+(y[combi[0]]-y[combi[dataLength-1]])**2)
        # それぞれの距離を足し合わせる
        totalDistance=disLast
        for dis in disz:
            totalDistance+=dis
        #最小距離が初期値or今のトータル距離より大きい時、データを上書き
        if minimumDis == 0 or minimumDis > totalDistance:
            minimumCombi=combi
            minimumDis=totalDistance
    return minimumCombi

def printResult(minimumCombi):  
    print("index")
    for number in minimumCombi:
        print(number)

x=[]
y=[]
csvToXyList()
indexs=makeIndexs()
minimumCombi=SearchMinimumDist(indexs)
printResult(minimumCombi)

    



    

