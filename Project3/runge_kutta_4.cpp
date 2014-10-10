#include<stdio.h>
//#include<conio.h>

float func(float x,float y);
void runge_kutta_4()
{
    float sx,ex,yp,h,i,yp1,s,s1,s2,s3,s4;
    
    // Variabelinput, må endres
    //clrscr();
    //printf("Enter the value of Y(1) ::");
    //scanf("%f",&yp);
    //printf("Enter the value of Start x ::");
    //scanf("%f",&sx);
    //printf("Enter the value of End x ::");
    //scanf("%f",&ex);
    //printf("Enter the value of h ::");
    //scanf("%f",&h);
    
    for(int i=sx;i<=ex+h;i=i+h)
    {
        s1 = func(i,yp);
        s2 = func(i+(h/2),yp+((h/2)*s1));
        s3 = func(i+(h/2),yp+((h/2)*s2));
        s4 = func(i+h,yp+(h*s3));
        s = (s1+(2*s2)+(2*s3)+s4)/6;
        yp1 = yp + (h*s);
        if((i*10)==(int)(ex*10))
            printf("\nx=%.4f\ty=%.4f",i,yp);
        else
            printf("\nx=%.4f\ty=%.4f\ty1=%.4f\th=%.4f\ts=%.4f",i,yp,yp1,h,s);
        yp = yp1;
    }
    getch();
}

float func(float x,float y)
{
    return x*y;
}

