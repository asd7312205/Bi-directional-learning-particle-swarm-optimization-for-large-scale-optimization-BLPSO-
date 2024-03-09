#include <iostream>
#include "CEC2013/Header.h"
#include <random>
#include <ctime>
#include <algorithm>
using namespace std;


mt19937 rng;
struct PARTICLE{
    double *position;
    double *velocity;
    double fit; //fitness
    int hIndex; //serial number of the bucket
    double dotProduct; //coordinates in one dimension
    int nh; //number of individuals in the same bucket
    double  excellentRatio;//convergence learning exempalr selection ratio
    double  sparseRatio;//diversity learning exempalr selection ratio
    int sortIndex;//ranking of individuals by fitnees
    double* pbestPosition;
    double pbestFit;

};

//fitness comparison function
bool cmp( PARTICLE i, PARTICLE j)
{
    return i.fit<j.fit;
}


//diversity comparison function
bool cmp1( PARTICLE i, PARTICLE j)
{
    return i.nh<j.nh;
}




Benchmarks* generateFuncObj(int funcID){
    Benchmarks *fp;
    // run each of specified function in "configure.ini"
    if (funcID==1){
        fp = new F1();
    }else if (funcID==2){
        fp = new F2();
    }else if (funcID==3){
        fp = new F3();
    }else if (funcID==4){
        fp = new F4();
    }else if (funcID==5){
        fp = new F5();
    }else if (funcID==6){
        fp = new F6();
    }else if (funcID==7){
        fp = new F7();
    }else if (funcID==8){
        fp = new F8();
    }else if (funcID==9){
        fp = new F9();
    }else if (funcID==10){
        fp = new F10();
    }else if (funcID==11){
        fp = new F11();
    }else if (funcID==12){
        fp = new F12();
    }else if (funcID==13){
        fp = new F13();
    }else if (funcID==14){
        fp = new F14();
    }else if (funcID==15){
        fp = new F15();
    }else{
        cerr<<"Fail to locate Specified Function Index"<<endl;
        exit(-1);
    }
    return fp;
}


void BLPSO(int fes, int initParticleSize, int dim)
{
    Benchmarks* fp=NULL; 
    double fai1,fai2;
    int runtimes;//FEs already used
    double w;
    int funToRun[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    int funNum = 15;
    int particleSize;
	
    ofstream write;
    double bestSize;//the size of convergence learning exempalrs
    double aerfa;//Maximum 2, minimum 0

    

        for (int cnum = 0;cnum < funNum; ++cnum) {
			if(cnum == 12 || cnum == 13)
            {
                dim=905;
            }
            else
            {
                dim=1000;
            }

            particleSize=initParticleSize;
            string name1="function";
            string name2;
            name2= to_string(cnum+1);
            name1+=name2;
            name2=".txt";
            name1+=name2;
            write.open(name1);

            runtimes=0;
            fp = generateFuncObj(funToRun[cnum]);
            double xMax=fp->getMaxX();
            double xmin=fp->getMinX();
            fp->setDimension(dim);


            rng.seed(random_device()());
            uniform_real_distribution<double> rd(xmin,xMax);
            uniform_real_distribution<double> rdv(xmin*0.2,xMax*0.2);
            uniform_real_distribution<double> rdr(0,1);
            uniform_real_distribution<double> rw(0,1);


            PARTICLE* pe=new PARTICLE [particleSize];
            for (int i = 0; i < particleSize; ++i) {
                pe[i].position=new double [dim];
                pe[i].velocity=new double [dim];
                pe[i].excellentRatio=0.1;
                pe[i].sparseRatio=0.1;
                pe[i].pbestPosition=new double [dim];

                for (int j = 0; j < dim; ++j) {
                    pe[i].position[j]=rd(rng);
                    pe[i].velocity[j]=rdv(rng);
                    pe[i].pbestPosition[j]=pe[i].position[j];

                }
                pe[i].fit=fp->compute(pe[i].position);
                runtimes++;
                pe[i].pbestFit=pe[i].fit;
            }


            double* ro=new double [dim];
            double hmax,hmin; //maximum and minimum values of coordinate positions
            double r;
            double nb=0.1*particleSize;

            int* countN=new int [(int)nb];
            int tr1,tr2,lshIndex;

            sort(pe,pe+particleSize,cmp);//sort by fitness
			
            bool* flags=new bool[6];
            for (int i = 0; i < 6; ++i) {
                flags[i]= true;
            }
            int* updateIndex=new int[particleSize];
            int countUpdateNumber;
            write<<runtimes<<" "<<pe[0].fit<<endl;
            int tempSize,tempiindex;

            while(runtimes < fes)
            {
                //LSH
                for (int i = 0; i < nb; ++i) {
                    countN[i]=0;
                }
                for (int i = 0; i < dim; ++i) {//randomized D-dimensional vectors
                    ro[i]=rd(rng);
                }
                for (int i = 0; i < particleSize; ++i) {
                    pe[i].dotProduct=0;
                    for (int j = 0; j < dim; ++j) {
                        pe[i].dotProduct+=pe[i].pbestPosition[j]*ro[j];
                    }
                }
                hmax=pe[0].dotProduct;
                hmin=hmax;
                for (int i = 0; i < particleSize; ++i) {
                    if(pe[i].dotProduct > hmax)
                    {
                        hmax=pe[i].dotProduct;
                    }
                    else if(pe[i].dotProduct < hmin)
                    {
                        hmin=pe[i].dotProduct;
                    }
                }

                r=(hmax-hmin)/nb;
                uniform_real_distribution<double> rb(0,r);
                for (int i = 0; i < particleSize; ++i) {
                    pe[i].hIndex= floor(((pe[i].dotProduct-hmin)+rb(rng))/r);
                    if(pe[i].hIndex < 0)
                    {
                        pe[i].hIndex=0;
                    }
                    if(pe[i].hIndex > nb-1)//forced division of the transgressions into the last bucket
                    {
                        pe[i].hIndex=nb-1;
                    }
                    countN[pe[i].hIndex]++;


                }
                for (int i = 0; i < particleSize; ++i) {
                    pe[i].nh=countN[pe[i].hIndex];

                }

                

                //record the individuals to be updated
                countUpdateNumber=0;
                for (int i = 0; i < particleSize; ++i) {
                    if((rdr(rng) < pow(((double)(i)/(double)(particleSize-1)),2)))
                    {
                        updateIndex[countUpdateNumber]=i;
                        countUpdateNumber++;
                    }
                }

                for (int i = countUpdateNumber-1; i > -1; --i) {
                    if(pe[updateIndex[i]].nh == 1)//individuals with the best diversity are retained
                    {
                        continue;
                    }

                    aerfa=(double)(updateIndex[i])/(double)(particleSize);
                    pe[updateIndex[i]].excellentRatio=1.0-aerfa;
                    bestSize= ceil(pe[updateIndex[i]].excellentRatio*updateIndex[i]);
                    fai1=1.0;
                    fai2=0.5-0.3*(double)(runtimes)/(double)(fes);
                    for (int j = 0; j < dim; ++j) {
                        tr1=rand()%(int)bestSize;//the index of the convergence learning exempalr 
                        w=rw(rng);
                        if(pe[updateIndex[i]].nh > 0.1*particleSize)
                        {
                            do {
                                tr2=rand()%particleSize;//the index of the diversity learning exempalr 
                            }while( (pe[tr2].nh > pe[updateIndex[i]].nh) );
                            pe[updateIndex[i]].velocity[j]=w*pe[updateIndex[i]].velocity[j]+fai1*rdr(rng)*(pe[tr1].position[j]-pe[updateIndex[i]].position[j])+fai2*rdr(rng)*(pe[tr2].pbestPosition[j]-pe[updateIndex[i]].position[j]);
                        }
                        else
                        {
                            tr2=rand()%(updateIndex[i]);
                            pe[updateIndex[i]].velocity[j]=w*pe[updateIndex[i]].velocity[j]+fai1*rdr(rng)*(pe[tr1].position[j]-pe[updateIndex[i]].position[j])+fai2*rdr(rng)*(pe[tr2].position[j]-pe[updateIndex[i]].position[j]);
                        }

                        pe[updateIndex[i]].position[j]=pe[updateIndex[i]].position[j]+ pe[updateIndex[i]].velocity[j];
                        if(pe[updateIndex[i]].position[j] > xMax)
                        {
                            pe[updateIndex[i]].position[j]=xMax;
                        }
                        else if(pe[updateIndex[i]].position[j] < xmin)
                        {
                            pe[updateIndex[i]].position[j]=xmin;
                        }
                    }
                    pe[updateIndex[i]].fit=fp->compute(pe[updateIndex[i]].position);
                    runtimes++;

                }

                //update the pbest
                for (int i = 0; i < countUpdateNumber; ++i) {
                    if(pe[updateIndex[i]].fit < pe[updateIndex[i]].pbestFit)
                    {
                        for (int j = 0; j < dim; ++j) {
                            pe[updateIndex[i]].pbestPosition[j]=pe[updateIndex[i]].position[j];
                        }
                        pe[updateIndex[i]].pbestFit=pe[updateIndex[i]].fit;
                    }
                }

                sort(pe,pe+particleSize,cmp);
                tempSize=800-round(600*(runtimes)/(fes));

                while(particleSize > tempSize)//random discard of redundant individuals
                {
                    tempiindex=(rand()%(particleSize-1))+1;
                    swap(pe[tempiindex],pe[particleSize-1]);
                    particleSize--;
                }
                //document the optimization results
                for (int i = 0; i < 6; ++i) {
                    if(runtimes > (i+1)*(fes/6) )
                    {
                        if(flags[i])
                        {
                            write<<runtimes<<" "<<pe[0].fit<<endl;
                            flags[i]= false;
                            break;
                        }
                        else
                        {
                            continue;
                        }
                    }
                    else
                    {
                        break;
                    }
                }

            }



            write<<endl;
            write<<"function"<<cnum+1<<": "<<pe[0].fit<<endl;
            write.close();
            delete fp;
            delete [] flags;
            delete [] countN;
            delete [] updateIndex;
            for (int i = 0; i < particleSize; ++i) {
                delete [] pe[i].position;
                delete [] pe[i].velocity;
                delete [] pe[i].pbestPosition;
            }

            delete [] pe;
        }

}

int main() {

    srand((unsigned )time(NULL));
    int dim = 1000;//dimension
    int particleSize=800;//the population size
    int fesMax=3000000;//the maximum of FEs
    
    BLPSO(fesMax,particleSize,dim);

    return 0;
}

