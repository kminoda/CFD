#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double re = 70.0;                                         //Reynolds
const double cfl = 0.2;                                         //CFL Number
const double omegap = 1.0;                                     //Blasius
const int maxitp = 100;
const double errorp = 1.0e-4;
const int nbegin = 0;
const int nlast = 5000;
// mesh
const int mx = 401;
const int my = 201;

// 正方形
const int i_1 = 96;
const int i_2 = 106;
const int j_1 = 91;
const int j_2 = 111;

const int icent = (i_1 + i_2)/2;
const int jcent = (j_1 + j_2)/2;


// 正方形の一辺の長さを元に計算格子の間隔を設定
const double dx = 1.0/(i_2-i_1);
const double dy = 2.0/(j_2-j_1);

// 時間間隔を計算
const double dt = cfl*min(dx,dy);

double resp,itrp;
double now_time = 0.0;

double x[mx + 2][my + 2], y[mx + 2][my + 2];
double u[mx + 2][my + 2], v[mx + 2][my + 2];
double p[mx + 2][my + 2];

double cp1,cp2;

/* print inti */
int print_init()
{
    printf("==================================================\n");
    printf("Incompressible Navier-Stokes 2D Flow Solver\n");
    printf("Flow around Rectangular Cylinder\n");
    printf("Programmed by K. Minoda (U-Tokyo AeroAstroEng. B3)\n");
    printf("==================================================\n\n");
    
    printf("*** Comp. Conditions\n");
    printf("CFL = %f\n",cfl);
    printf("dt = %f\n",dt);
    printf("%d Time Steps to go...\n\n",nlast);
    printf(">>2D Incompressible Flow Solver\n");
    printf("RE = %f\n",re);
    printf("No. of Grid Points: %d,%d\n",mx,my);
    return 0;
}

// 計算格子の設定
void setgrd()
{
    for (int i = 1; i < mx+1 ; i++)
    {
        for (int j = 1; j < my+1 ; j++)
        {
            x[i][j] = dx * (i-icent);
            y[i][j] = dy * (j-jcent);
        }
    }
}

// 初期条件の設定
void intcnd()
{
    for (int i = 1; i < mx ; i++)
    {
        for (int j = 1; j < my ; j++)
        {
            u[i][j] = 1.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
        }
    }
}

// 圧力の境界条件の設定
void bcforp()
{
    for (int j = 1; j <= my ; j++)
    {
        p[1][j] = 0.0;  // inflow condition ( i=1 )
        p[mx][j] = 0.0; // downstream condition ( i=mx )
    }
    for (int i = 1; i <= mx ; i++)
    {
        p[i][1] = 0.0;  // bottom condition ( j=1 )
        p[i][my] = 0.0; // bottom contition ( j=my )
    }
    // wall condition
    p[i_1][j_1] = p[i_1 - 1][j_1 - 1];
    p[i_1][j_2] = p[i_1 - 1][j_2 + 1];
    p[i_2][j_1] = p[i_2 + 1][j_1 - 1];
    p[i_2][j_2] = p[i_2 + 1][j_2 + 1];
    for (int j = j_1+1; j <= j_2-1 ; j++)
    {
        p[i_1][j] = p[i_1-1][j];
        p[i_2][j] = p[i_2+1][j];
    }
    for (int i = i_1+1; i <= i_2-1 ; i++)
    {
        p[i][j_1] = p[i][j_1-1];
        p[i][j_2] = p[i][j_2+1];
    }
}


void poiseq()
{
    double rhs[mx + 1][my + 1];
    double dp;
    double ux,uy,vx,vy;
    double res;
    int itr;
    for (int i = 2; i <= mx-1 ; i++)
    {
        for (int j = 2; j <= my-1 ; j++)
        {
            if( (i>i_1 && i<i_2) && (j>=j_1 && j<=j_2) ){continue;}
            ux = (u[i+1][j]-u[i-1][j])/(2.0*dx);
            uy = (u[i][j+1]-u[i][j-1])/(2.0*dy);
            vx = (v[i+1][j]-v[i-1][j])/(2.0*dx);
            vy = (v[i][j+1]-v[i][j-1])/(2.0*dy);
            rhs[i][j] = (ux+vy)/dt - (ux*ux + 2.0*uy*vx + vy*vy);
        }
    }
    
    for (itr = 1; itr <= maxitp; itr++)
    {
        double res = 0.0;   // residual（残差）
        for (int i = 2; i <= mx-1 ; i++)
        {
            for (int j = 2; j <= my-1 ; j++)
            {
            if( (i>i_1 && i<i_2) && (j>=j_1 && j<=j_2) ){continue;}
                dp = (p[i+1][j]+p[i-1][j])/(dx*dx) + (p[i][j+1]+p[i][j-1])/(dy*dy) -rhs[i][j];
                dp = dp/(2.0/(dx*dx) + 2.0/(dy*dy))-p[i][j];
                res = res + dp*dp;
                p[i][j] = p[i][j] + omegap*dp;
            }
        }
        bcforp();
        res = pow(res/double(mx*my),0.5);
        
        //printf("res = %f ( at itr = %d )\n",res,itr);
        
        
        if(res < errorp)
        {
            resp = res;
            itrp = itr; // itr;
            break;
        }
    }
}

void veloeq()
{
    double urhs[mx + 1][my + 1];
    double vrhs[mx + 1][my + 1];
    for (int i = 2; i <= mx-1 ; i++)
    {
        for (int j = 2; j <= my-1 ; j++)
        {
            if( (i>i_1 && i<i_2) && (j>=j_1 && j<=j_2) ){continue;}
            urhs[i][j] = -(p[i+1][j] - p[i-1][j])/(2.0*dx);
            vrhs[i][j] = -(p[i][j+1] - p[i][j-1])/(2.0*dy);
        }
    }
    for (int i = 2; i <= mx-1 ; i++)
    {
        for (int j = 2; j <= my-1 ; j++)
        {
            if( (i>i_1 && i<i_2) && (j>=j_1 && j<=j_2) ){continue;}
            urhs[i][j] = urhs[i][j] + (u[i+1][j]-2.0*u[i][j]+u[i-1][j])/(re*dx*dx) +(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/(re*dy*dy);
            vrhs[i][j] = vrhs[i][j] + (v[i+1][j]-2.0*v[i][j]+v[i-1][j])/(re*dx*dx) + (v[i][j+1]-2.0*v[i][j]+v[i][j-1])/(re*dy*dy);
        }
    }
    for (int j = j_1+1; j <= j_2-1 ; j++)
    {
        u[i_1+1][j] = 2.0*u[i_1][j]-u[i_1-1][j];
        u[i_2-1][j] = 2.0*u[i_2][j]-u[i_2+1][j];
        v[i_1+1][j] = 2.0*v[i_1][j]-v[i_1-1][j];
        v[i_2-1][j] = 2.0*v[i_2][j]-v[i_2+1][j];
    }
    for (int i = 2; i <= mx-1 ; i++)
    {
        for (int j = 2; j <= my-1 ; j++)
        {
            if( (i>i_1 && i<i_2) && (j>=j_1 && j<=j_2) ){continue;}
            urhs[i][j] = urhs[i][j] - u[i][j]*(-u[i+2][j] + 8.0*(u[i+1][j]-u[i-1][j])+u[i-2][j])/(12.0*dx) - abs(u[i][j])*(u[i+2][j]-4.0*u[i+1][j]+6.0*u[i][j]-4.0*u[i-1][j]+u[i-2][j])/(4.0*dx);
            vrhs[i][j] = vrhs[i][j] - u[i][j]*(-v[i+2][j] + 8.0*(v[i+1][j]-v[i-1][j])+v[i-2][j])/(12.0*dx) - abs(u[i][j])*(v[i+2][j]-4.0*v[i+1][j]+6.0*v[i][j]-4.0*v[i-1][j]+v[i-2][j])/(4.0*dx);
        }
    }
    for (int i = i_1+1; i <= i_2-1 ; i++)
    {
        u[i][j_1+1] = 2.0*u[i][j_1]-u[i][j_1-1];
        u[i][j_2-1] = 2.0*u[i][j_2]-u[i][j_2+1];
        v[i][j_1+1] = 2.0*v[i][j_1]-v[i][j_1-1];
        v[i][j_2-1] = 2.0*v[i][j_2]-v[i][j_2+1];
    }
    for (int i = 2; i <= mx-1 ; i++)
    {
        for (int j = 2; j <= my-1 ; j++)
        {
            if( (i>i_1 && i<i_2) && (j>=j_1 && j<=j_2) ){continue;}
            urhs[i][j] = urhs[i][j] - v[i][j]*(-u[i][j+2] + 8.0*(u[i][j+1]-u[i][j-1])+u[i][j-2])/(12.0*dy) - abs(v[i][j])*(u[i][j+2]-4.0*u[i][j+1]+6.0*u[i][j]-4.0*u[i][j-1]+u[i][j-2])/(4.0*dy);
            vrhs[i][j] = vrhs[i][j] - v[i][j]*(-v[i][j+2] + 8.0*(v[i][j+1]-v[i][j-1])+v[i][j-2])/(12.0*dy) - abs(v[i][j])*(v[i][j+2]-4.0*v[i][j+1]+6.0*v[i][j]-4.0*v[i][j-1]+v[i][j-2])/(4.0*dy);
        }
    }
    for (int i = 2; i <= mx-1 ; i++)
    {
        for (int j = 2; j <= my-1 ; j++)
        {
            if( (i>i_1 && i<i_2) && (j>=j_1 && j<=j_2) ){continue;}
            u[i][j] = u[i][j]+dt*urhs[i][j];
            v[i][j] = v[i][j]+dt*vrhs[i][j];
        }
    }
}

void bcforv()
{
    for (int j = 1; j <= my ; j++)
    {
        u[1][j] = 1.0;
        v[1][j] = 0.0;
        u[0][j] = 0.0;
        v[0][j] = 0.0;
        u[mx][j] = 2.0*u[mx-1][j] - u[mx-2][j];
        v[mx][j] = 2.0*v[mx-1][j] - v[mx-2][j];
        u[mx+1][j] = 2.0*u[mx][j] - u[mx-1][j];
        v[mx+1][j] = 2.0*v[mx][j] - v[mx-1][j];
    }
    for (int i = 1; i <= mx ; i++)
    {
        u[i][1] = 2.0*u[i][2] - u[i][3];
        v[i][1] = 2.0*v[i][2] - v[i][3];
        u[i][0] = 2.0*u[i][1] - u[i][2];
        v[i][0] = 2.0*v[i][1] - v[i][2];
        u[i][my] = 2.0*u[i][my-1] - u[i][my-2];
        v[i][my] = 2.0*v[i][my-1] - v[i][my-2];
        u[i][my+1] = 2.0*u[i][my] - u[i][my-1];
        v[i][my+1] = 2.0*v[i][my] - v[i][my-1];
    }
    // wall condition
    for (int i = i_1; i <= i_2 ; i++)
    {
        for (int j = j_1; j <= j_2 ; j++)
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
}

void slvflw()
{
    printf("now we are at slvflw\n");
    double nstep;
    int nlp =10;
    int nlp_csv = 2500;
    printf("%f\n",dt);
    intcnd();
    bcforv();
    bcforp();
    printf("Step / Res(p) at itr. / CD / CL / Cp1 / Cp2\n");
    
    /*
    std::fstream ResFile500("result_t500.csv", std::ios::out);
    std::fstream ResFile500("result_t1000.csv", std::ios::out);
    std::fstream ResFile500("result_t1500.csv", std::ios::out);
    std::fstream ResFile500("result_t2000.csv", std::ios::out);
    std::fstream ResFile500("result_t2500.csv", std::ios::out);
    std::fstream ResFile500("result_t3000.csv", std::ios::out);
    std::fstream ResFile500("result_t3500.csv", std::ios::out);
    std::fstream ResFile500("result_t4000.csv", std::ios::out);
    std::fstream ResFile500("result_t4500.csv", std::ios::out);
    //ResFile << "time" << "," << "x" << "," << "y" << "," << "u" << "," << "v" << "," << "p" <<  "," << "cd" <<  "," << "cl" << "," << "cp1" << "," << "cp2" << std::endl;
    ResFile500 << "x"<<","<<"y"<<","<<"p" << std::endl;
    ResFile1000 << "x"<<","<<"y"<<","<<"p" << std::endl;
    ResFile1500 << "x"<<","<<"y"<<","<<"p" << std::endl;
    ResFile2000 << "x"<<","<<"y"<<","<<"p" << std::endl;
    ResFile2500 << "x"<<","<<"y"<<","<<"p" << std::endl;
    ResFile3000 << "x"<<","<<"y"<<","<<"p" << std::endl;
    ResFile3500 << "x"<<","<<"y"<<","<<"p" << std::endl;
    ResFile4000 << "x"<<","<<"y"<<","<<"p" << std::endl;
    ResFile4500 << "x"<<","<<"y"<<","<<"p" << std::endl; */
    for (int n = 1; n <= nlast; n++)
    {
        nstep = n+nbegin;
        now_time = now_time + dt;
        poiseq();
        bcforp();
        veloeq();
        bcforv();
        double cd=0.0;
        double cpfore,cpback;
        for (int j = j_1; j <= j_2 - 1 ; j++)
        {
            cpfore = (2.0*p[i_1][j] + 2.0+p[i_1][j+1])/2.0;
            cpback = (2.0*p[i_2][j] + 2.0+p[i_2][j+1])/2.0;
            cd = cd + (cpfore-cpback)*dy;
        }
        double cl=0.0;
        double cpbtm,cptop;
        for (int i = i_1; i <= i_2 - 1 ; i++)
        {
            cpbtm = (2.0*p[i][j_1] + 2.0*p[i+1][j_1])/2.0;
            cptop = (2.0*p[i][j_2] + 2.0*p[i+1][j_2])/2.0;
            cl = cl + (cpbtm-cptop)*dx;
        }
        cp1 = 2.0*p[i_2+i_2-i_1][j_1];
        cp2 = 2.0*p[i_2+i_2-i_1][j_2];

        
        if(n==500)
        {
            std::fstream ResFile500("result_t500.csv", std::ios::out);
            ResFile500 << "x"<<","<<"y"<<","<<"p" << std::endl;
            for (int i = 1; i <= mx ; i++)
            {
                for (int j = 1; j <= my ; j++)
                {
                    //ResFile <<  now_time << "," << x[i][j] << "," << y[i][j] << "," << u[i][j] << "," << v[i][j] << "," << p[i][j] <<  "," << cd <<  "," << cl << "," << cp1 << "," << cp2 << std::endl;
                    ResFile500 << x[i][j] << "," << y[i][j] << "," << p[i][j] << std::endl;
                }
                ResFile500 << "" << "," << "" << "," << "" << std::endl;
            }
        }

        if(n==1000)
        {
            std::fstream ResFile1000("result_t1000.csv", std::ios::out);
            ResFile1000 << "x"<<","<<"y"<<","<<"p" << std::endl;
            for (int i = 1; i <= mx ; i++)
            {
                for (int j = 1; j <= my ; j++)
                {
                    //ResFile <<  now_time << "," << x[i][j] << "," << y[i][j] << "," << u[i][j] << "," << v[i][j] << "," << p[i][j] <<  "," << cd <<  "," << cl << "," << cp1 << "," << cp2 << std::endl;
                    ResFile1000 << x[i][j] << "," << y[i][j] << "," << p[i][j] << std::endl;
                }
                ResFile1000 << "" << "," << "" << "," << "" << std::endl;
            }
        }
        if(n==1500)
        {
            std::fstream ResFile1500("result_t1500.csv", std::ios::out);
            ResFile1500 << "x"<<","<<"y"<<","<<"p" << std::endl;
            for (int i = 1; i <= mx ; i++)
            {
                for (int j = 1; j <= my ; j++)
                {
                    //ResFile <<  now_time << "," << x[i][j] << "," << y[i][j] << "," << u[i][j] << "," << v[i][j] << "," << p[i][j] <<  "," << cd <<  "," << cl << "," << cp1 << "," << cp2 << std::endl;
                    ResFile1500 << x[i][j] << "," << y[i][j] << "," << p[i][j] << std::endl;
                }
                ResFile1500 << "" << "," << "" << "," << "" << std::endl;
            }
        }
        if(n==2000)
        {
            std::fstream ResFile2000("result_t2000.csv", std::ios::out);
            ResFile2000 << "x"<<","<<"y"<<","<<"p" << std::endl;
            for (int i = 1; i <= mx ; i++)
            {
                for (int j = 1; j <= my ; j++)
                {
                    //ResFile <<  now_time << "," << x[i][j] << "," << y[i][j] << "," << u[i][j] << "," << v[i][j] << "," << p[i][j] <<  "," << cd <<  "," << cl << "," << cp1 << "," << cp2 << std::endl;
                    ResFile2000 << x[i][j] << "," << y[i][j] << "," << p[i][j] << std::endl;
                }
                ResFile2000 << "" << "," << "" << "," << "" << std::endl;
            }
        }
        if(n==2500)
        {
            std::fstream ResFile2500("result_t2500.csv", std::ios::out);
            ResFile2500 << "x"<<","<<"y"<<","<<"p" << std::endl;
            for (int i = 1; i <= mx ; i++)
            {
                for (int j = 1; j <= my ; j++)
                {
                    //ResFile <<  now_time << "," << x[i][j] << "," << y[i][j] << "," << u[i][j] << "," << v[i][j] << "," << p[i][j] <<  "," << cd <<  "," << cl << "," << cp1 << "," << cp2 << std::endl;
                    ResFile2500 << x[i][j] << "," << y[i][j] << "," << p[i][j] << std::endl;
                }
                ResFile2500 << "" << "," << "" << "," << "" << std::endl;
            }
        }
        if(n==3000)
        {
            std::fstream ResFile3000("result_t3000.csv", std::ios::out);
            ResFile3000 << "x"<<","<<"y"<<","<<"p" << std::endl;
            for (int i = 1; i <= mx ; i++)
            {
                for (int j = 1; j <= my ; j++)
                {
                    //ResFile <<  now_time << "," << x[i][j] << "," << y[i][j] << "," << u[i][j] << "," << v[i][j] << "," << p[i][j] <<  "," << cd <<  "," << cl << "," << cp1 << "," << cp2 << std::endl;
                    ResFile3000 << x[i][j] << "," << y[i][j] << "," << p[i][j] << std::endl;
                }
                ResFile3000 << "" << "," << "" << "," << "" << std::endl;
            }
        }
        if(n==3500)
        {
            std::fstream ResFile3500("result_t3500.csv", std::ios::out);
            ResFile3500 << "x"<<","<<"y"<<","<<"p" << std::endl;
            for (int i = 1; i <= mx ; i++)
            {
                for (int j = 1; j <= my ; j++)
                {
                    //ResFile <<  now_time << "," << x[i][j] << "," << y[i][j] << "," << u[i][j] << "," << v[i][j] << "," << p[i][j] <<  "," << cd <<  "," << cl << "," << cp1 << "," << cp2 << std::endl;
                    ResFile3500 << x[i][j] << "," << y[i][j] << "," << p[i][j] << std::endl;
                }
                ResFile3500 << "" << "," << "" << "," << "" << std::endl;
            }
        }
        if(n==4000)
        {
            std::fstream ResFile4000("result_t4000.csv", std::ios::out);
            ResFile4000 << "x"<<","<<"y"<<","<<"p" << std::endl;
            for (int i = 1; i <= mx ; i++)
            {
                for (int j = 1; j <= my ; j++)
                {
                    //ResFile <<  now_time << "," << x[i][j] << "," << y[i][j] << "," << u[i][j] << "," << v[i][j] << "," << p[i][j] <<  "," << cd <<  "," << cl << "," << cp1 << "," << cp2 << std::endl;
                    ResFile4000 << x[i][j] << "," << y[i][j] << "," << p[i][j] << std::endl;
                }
                ResFile4000 << "" << "," << "" << "," << "" << std::endl;
            }
        }
        if(n==4500)
        {
            std::fstream ResFile4500("result_t4500.csv", std::ios::out);
            ResFile4500 << "x"<<","<<"y"<<","<<"p" << std::endl;
            for (int i = 1; i <= mx ; i++)
            {
                for (int j = 1; j <= my ; j++)
                {
                    //ResFile <<  now_time << "," << x[i][j] << "," << y[i][j] << "," << u[i][j] << "," << v[i][j] << "," << p[i][j] <<  "," << cd <<  "," << cl << "," << cp1 << "," << cp2 << std::endl;
                    ResFile4500 << x[i][j] << "," << y[i][j] << "," << p[i][j] << std::endl;
                }
                ResFile4500 << "" << "," << "" << "," << "" << std::endl;
            }
        }

        if(!(n%nlp))
        {
            printf("%d / %d / %f / %f / %f / %f\n",int(nstep),int(itrp),cd,cl,cp1,cp2);
        }
    }
    for (int i = 1; i <= mx ; i++)
    {
        for (int j = 1; j <= my ; j++)
        {
            p[i][j] = 2.0*p[i][j];
        }
    }
}


int main()
{
    print_init();
    setgrd();
    slvflw();

	return 0;
}
