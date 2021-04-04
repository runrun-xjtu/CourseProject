#include <iostream>
#include <iomanip>
#include <math.h>
#include <ctime>
using namespace std;

/*可调公式与变量区域*/
#define f(x1,x2) x1*x1+x2;          //f(x1,x2) 表达式
#define g1(x1,x2) x1*x1+x2*x2-9;    //g1(x1,x2) 表达式
#define g2(x1,x2) x1+x2-1;          //g2(x1,x2) 表达式
const double x1Max = 3;				//x1 取值区间
const double x1Min = -3;
const double x2Max = 3;				//x2 取值区间
const double x2Min = -3;
const int popSize = 20;		        //*种群大小
const int len = 10;					//*染色体长度（控制精度）
const double Pc = 0.9;		        //*交叉概率
const double Pm = 0.1;		        //*变异概率
const int generation = 50;         //*迭代次数（终止条件）

typedef struct individual   //个体结构体
{
	bool chromo1[len];      //染色体 x1
	bool chromo2[len];      //染色体 x2
	double fit;             //适应度
	double prob;            //被选中概率
}individual;
struct individual pop[popSize];         //当代种群
struct individual nextPop[popSize];     //下一代种群

double totalFit;                   //总适应度
int gen = 1;                       //当前代数
double x1Range = x1Max - x1Min;    //x1 区间长度
double x2Range = x2Max - x2Min;    //x2 区间长度
double fxmin = 1000000;            //f(x)的最小值（所有个体中出现过的最小值）,初值 1000000
int nextPopNum = 0;                //下一代已产生个体的计数

void initPop(struct individual pop[]);     
void evolution();
void selection(struct individual[]);
void crossover(struct individual[]);
void mutataion(struct individual pop[]);

bool inLegality(individual);
void updateFmin(struct individual pop[]);
void getFitness(struct individual pop[]);
double getTotalFit(struct individual pop[]);
double getProb(individual);
double getValue(bool chromo[]);
void popPrintf(struct individual pop[]);

int main()
{
	double x1opti, x2opti, fxopti;
	double x1result[100], x2result[100], fresult[100];
	int n = 5;
    srand((unsigned)time(NULL));

	for (int i = 0; i < n; i++)
	{
		gen = 1; fxmin = 1000000; nextPopNum = 0;
		initPop(pop);
		while (gen < generation)
			evolution();
		popPrintf(pop);
		x1opti = x1Range * (getValue(pop[0].chromo1) / (pow(2, len) - 1)) + x1Min;
		x2opti = x2Range * (getValue(pop[0].chromo2) / (pow(2, len) - 1)) + x2Min;
		fxopti = f(x1opti, x2opti);
		cout << "最优解为：f(" << x1opti << "," << x2opti << ")=" << fxopti << endl;
		x1result[i] = x1opti;
		x2result[i] = x2opti;
		fresult[i] = fxopti;
	}
	for (int i = 0; i < n; i++)
	{
		cout <<"试验" << i+1 <<" 最优解为：f(" << x1result[i] << "," << x2result[i] << ")=     " << fresult[i] << endl;
	}
	return 0;
}

/*生成随机初始种群*/
void initPop(struct individual pop[])
{
	double g1value, g2value;
	double x1, x2;
	int	x1Num,x2Num;
	int i = 0;

	while (i < popSize)
	{
		x1 = x1Range * (rand() / double(RAND_MAX)) + x1Min;      //生成 -3 到 3 的随机数赋给x1
		x2 = x2Range * (rand() / double(RAND_MAX)) + x2Min;      //生成 -3 到 3 的随机数赋给x2
		g1value = g1(x1, x2);
		g2value = g2(x1, x2);
		if (g1value <= 0 && g2value <= 0)                                  //如果随机坐标（x1,x2）落入可行域，则取用
		{
			x1Num = floor((x1 - x1Min) / x1Range * (pow(2, len) - 1));      //获得 x1 对应的十进制排序值（0 - 2^len-1区间）
			x2Num = floor((x2 - x1Min) / x1Range * (pow(2, len) - 1));      //获得 x2 对应的十进制排序值（0 - 2^len-1区间）
			for (int j = 0; j < len; j++)
			{
				pop[i].chromo1[j] = floor(x1Num / (pow(2, len - j - 1)));   //x1 染色体赋值
				x1Num = x1Num % (int)pow(2, len - j - 1);
				pop[i].chromo2[j] = floor(x2Num / (pow(2, len - j - 1)));   //x2 染色体赋值
				x2Num = x2Num % (int)pow(2, len - j - 1);
			}
			i++;
		}
	}
	updateFmin(pop);                          //更新f(x1,x2)的最小值
	getFitness(pop);                          //更新个体适应度
	totalFit = getTotalFit(pop);              //更新总适应度
	for (int i = 0; i < popSize; i++)         
		pop[i].prob = getProb(pop[i]);        //更新选择概率
	cout << "第1代种群：" << endl;	          //打印初始种群
	popPrintf(pop);
}

/*进化到下一代*/
void evolution()
{
	nextPopNum = 0;     //下一代个体计数
	gen++;
	cout << "第" << gen << "代种群：" << endl;

	while (nextPopNum < popSize)
	{
		selection(pop);	                         //复制1个个体到下一代
		if (nextPopNum == popSize) break;
		crossover(pop);                          //生成1/2个杂交个体到下一代
		if (nextPopNum == popSize) break;
		mutataion(pop);                          //变异1个个体到下一代 
		if (nextPopNum == popSize) break;
	}
	updateFmin(nextPop);                             //更新f(x1,x2)的最小值
	getFitness(nextPop);                             //更新个体适应度
	totalFit = getTotalFit(nextPop);             //更新适应度与选择频率
	for (int i = 0; i < popSize; i++)
		nextPop[i].prob = getProb(nextPop[i]);
	popPrintf(nextPop);

	for (int i = 0; i < popSize; i++)            //更新种群存储空间
		pop[i] = nextPop[i];
}

/*选择以复制（轮盘式）*/
void selection(struct individual pop[])
{
	double lineProb[popSize];
	double r = rand() / double(RAND_MAX);           //获取0~1之间随机数
	cout << "  选择一次" << endl;
	lineProb[0] = pop[0].prob;
	for (int i = 1; i < popSize; i++)               //累计概率，使种群内个体的概率线性化
	{
		lineProb[i] = lineProb[i - 1] + pop[i].prob;
	}
	for (int j = 0; j < popSize; j++)               //查询随机数对应的轮盘上的个体
	{
		if (r <= lineProb[j])
		{
			nextPop[nextPopNum] = pop[j];
			nextPopNum++;
			//cout << nextPopNum << endl;
			break;
		}
	}
}

/*杂交*/
void crossover(struct individual pop[])
{
	double r1;
	int r2, r3, r4;
	bool a, b;
	r1 = rand() / double(RAND_MAX);//获取0~1之间随机数
	if (r1 <= Pc)
	{
		r2 = (rand() % (popSize));
		r3 = (rand() % (popSize));
		r4 = (rand() % (len));

		if (nextPopNum > popSize - 2)//当前新种群只剩一个位置
		{
			cout << "  杂交一次，产生一个新个体" << endl;
			nextPop[nextPopNum] = pop[r2];
			for (int i = r4; i < len; i++)
			{
				nextPop[nextPopNum].chromo1[i] = pop[r3].chromo1[i];
				nextPop[nextPopNum].chromo2[i] = pop[r3].chromo2[i];
			}
			if (inLegality(nextPop[nextPopNum])) //判断新个体值是否合法
				nextPopNum++;
		}
		else //新种群中剩余空位大于等于2
		{
			cout << "  杂交一次，产生二个新个体" << endl;
			nextPop[nextPopNum] = pop[r2];
			nextPopNum++;
			nextPop[nextPopNum] = pop[r3];

			//两个个体进行杂交
			for (int i = r4; i < len; i++)
			{
				a = nextPop[nextPopNum - 1].chromo1[i];
				nextPop[nextPopNum - 1].chromo1[i] = nextPop[nextPopNum].chromo1[i];
				nextPop[nextPopNum].chromo1[i] = a;

				b = nextPop[nextPopNum - 1].chromo2[i];
				nextPop[nextPopNum - 1].chromo2[i] = nextPop[nextPopNum].chromo2[i];
				nextPop[nextPopNum].chromo2[i] = b;
			}
			//判断新个体值是否合法
			if (!inLegality(nextPop[nextPopNum - 1]))
			{
				if (!inLegality(nextPop[nextPopNum]))
					nextPopNum--;
				else
					nextPop[nextPopNum - 1] = nextPop[nextPopNum];
			}
			else
			{
				if (inLegality(nextPop[nextPopNum]))
					nextPopNum++;
			}
		}
	}
	else
		cout << "  此次未进行杂交" << endl;
}

/*变异*/
void mutataion(struct individual pop[])
{
	double m0, m1;
	int m2, m3;
	m0 = rand() / double(RAND_MAX);    //获取0~1之间随机数
	m1 = (rand() % (3));               //获取0、1之间随机数
	if (m0 < Pm)
	{
		cout << "  变异一次，产生一个新个体" << endl;
		m2 = (rand() % (popSize));
		m3 = (rand() % (len));
		nextPop[nextPopNum] = pop[m2];

		if (m1 == 0)
		{
			if (nextPop[nextPopNum].chromo1[m3] == 1)   //取反
				nextPop[nextPopNum].chromo1[m3] = 0;
			else
				nextPop[nextPopNum].chromo1[m3] = 1;
		}
		else if(m1 == 1)
		{
			if (nextPop[nextPopNum].chromo2[m3] == 1)   //取反
				nextPop[nextPopNum].chromo2[m3] = 0;
			else
				nextPop[nextPopNum].chromo2[m3] = 1;
		}
		else
		{
			if (nextPop[nextPopNum].chromo1[m3] == 1)   //取反
				nextPop[nextPopNum].chromo1[m3] = 0;
			else
				nextPop[nextPopNum].chromo1[m3] = 1;
			if (nextPop[nextPopNum].chromo2[m3] == 1)   //取反
				nextPop[nextPopNum].chromo2[m3] = 0;
			else
				nextPop[nextPopNum].chromo2[m3] = 1;
		}
		//判断新个体值是否合法
		if (inLegality(nextPop[nextPopNum]))
			nextPopNum++;
	}
	else
		cout << "  此次未进行变异" << endl;
}


/*个体是否在可行域内判断*/
bool inLegality(individual in)
{
	double x1, x2;
	double g1value, g2value;
	x1 = x1Range * (getValue(in.chromo1) / (pow(2, len) - 1)) + x1Min;
	x2 = x2Range * (getValue(in.chromo2) / (pow(2, len) - 1)) + x2Min;
	g1value = g1(x1, x2);
	g2value = g2(x1, x2);
	if( g1value <= 0 && g2value <= 0)
		return 1;
	else
		return 0;
}

/*更新f(x1,x2)最小值*/
void updateFmin(struct individual pop[])
{
	double x1, x2;
	double fvalue;
	for (int i = 0; i < popSize; i++)
	{
		x1 = x1Range * (getValue(pop[i].chromo1) / (pow(2, len) - 1)) + x1Min;
		x2 = x2Range * (getValue(pop[i].chromo2) / (pow(2, len) - 1)) + x2Min;
		fvalue = f(x1, x2);
		if (fvalue < fxmin)
		    fxmin = fvalue;
	}
}

/*获取个体适应度*/
void getFitness(struct individual pop[])
{
	double FitnessMax = 0;   
	double x1, x2, fvalue;
	double fit = 0;
	int EncourageToSurvive = 3 * popSize, numOfMin = 0;    //鼓励存活系数、达到“近收敛值”的个体数量（鼓励收敛时有利变异的个体存活）
	for (int i = 0; i < popSize; i++)
	{
		x1 = x1Range * (getValue(pop[i].chromo1) / (pow(2, len) - 1)) + x1Min;
		x2 = x2Range * (getValue(pop[i].chromo2) / (pow(2, len) - 1)) + x2Min;
		fvalue = f(x1, x2);
		if (fvalue == fxmin)
			pop[i].fit = 0;                        //当fvalue == fxmin时，适应度无法计算，暂时取为0
		else
			pop[i].fit = 1 / (fvalue - fxmin);     //更新个体适应度
		if (pop[i].fit > FitnessMax)
			FitnessMax = pop[i].fit;               //更新个体适应度最大值
	}
	for (int i = 0; i < popSize; i++)
	{
		for (int j = 0; j < 3; j++)               //判断是否达到近收敛值
		{
			if ((pop[i].chromo1 == pop[(rand() % (popSize))].chromo1) && (pop[i].chromo2 == pop[(rand() % (popSize))].chromo2))
			{
				numOfMin++;
				break;
			}
		}
		if (pop[i].fit == 0)
		{
			if ((numOfMin / popSize) <= 0.8)
			pop[i].fit = 5 * FitnessMax;                      //如果某个体适应度为0（当fvalue == fxmin时），将该值设置为该种群最大适应度的5倍
			else
			pop[i].fit = EncourageToSurvive * FitnessMax;     //如果某个体适应度为0且处于近收敛状态，将该值设置为该种群最大适应度的“鼓励存活系数”倍
		}		           
	}
}

/*获取总适应度*/
double getTotalFit(struct individual pop[])
{
	double totalFit = 0;
	for (int i = 0; i < popSize; i++)
		totalFit += pop[i].fit;
	return totalFit;
}

/*获取个体选择概率*/
double getProb(individual in)
{
	double prob = 0;
	prob = in.fit / totalFit;
	return prob;
}

/*获取二进制编码对应十进制值*/
double getValue(bool chromo[len])
{
	double x = 0;
	for (int i = 0; i < len; i++)
	{
		x += chromo[i] * pow(2, (len - i - 1));
	}
	return x;
}

/*打印种群*/
void popPrintf(struct individual pop[])
{
	double x1[popSize], x2[popSize],fvalue;
	for (int i = 0; i < popSize; i++)
	{
		x1[i] = x1Range * (getValue(pop[i].chromo1) / (pow(2, len) - 1)) + x1Min;
		x2[i] = x2Range * (getValue(pop[i].chromo2) / (pow(2, len) - 1)) + x2Min;
		fvalue = f(x1[i], x2[i]);
		
		cout << setw(5) << left << i;
		cout << setw(12) << left << x1[i];
		cout << setw(12) << left << x2[i];
		cout << "f(x1,x2) = ";
		cout << setw(12) << left << fvalue << "   ";
		for (int j = 0; j < len; j++)
			cout << pop[i].chromo1[j];
		cout << " ";
		for (int j = 0; j < len; j++)
			cout << pop[i].chromo2[j];
		cout << endl;
	}
}
