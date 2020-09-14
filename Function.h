#include "Constants.h"
class Function {
private:
    float x, y, r;
    metodo f1;
    float c1[2], c2[2];
    int dx, dy;
    float m;
    float steps;
    float pn;
    float Xinc, Yinc;
public:
    GLvoid graphing_circle() {

        
        r = sqrt(pow(c1[0], 2) + pow(c1[1], 2));
        for (float i = 0; i < r; i += 0.05) {
            y = sqrt(pow(r, 2) - pow(i, 2));
            f1.Pinta_Pixel(i, y, 71, 178,227);
        }
        for (float i = 0; i < r; i += 0.05) {
            y = -sqrt(pow(r, 2) - pow(i, 2));
            f1.Pinta_Pixel(i, y, 71, 178, 227);
        }
        for (float i = 0; i < r; i += 0.05) {
            y = sqrt(pow(r, 2) - pow(i, 2));
            f1.Pinta_Pixel(-i, y, 71, 178, 227);
        }
        for (float i = 0; i < r; i += 0.05) {
            y = -sqrt(pow(r, 2) - pow(i, 2));
            f1.Pinta_Pixel(-i, y, 71, 178, 227);
        }
        
        
    }
    GLvoid graphing_elipse() {

        float r1 = sqrt(pow(c1[0], 2) + pow(c1[1],2));
        float r2 = sqrt(pow(c2[0], 2) + pow(c2[1], 2));
        
        for (float i = 0; i < r1; i += 0.05) {

            y = sqrt(pow(r2, 2) * (1 - (pow(i, 2) / pow(r1, 2))));
            f1.Pinta_Pixel(i, y, 1, 1, 1);
        }
        for (float i = 0; i < r1; i += 0.05) {

            y = -sqrt(pow(r2, 2) * (1 - (pow(i, 2) / pow(r1, 2))));
            f1.Pinta_Pixel(i, y, 1, 1, 1);
        }
        for (float i = 0; i < r1; i += 0.05) {

            y = sqrt(pow(r2, 2) * (1 - (pow(i, 2) / pow(r1, 2))));
            f1.Pinta_Pixel(-i, y, 1, 1, 1);
        }
        for (float i = 0; i < r1; i += 0.05) {

            y = -sqrt(pow(r2, 2) * (1 - (pow(i, 2) / pow(r1, 2))));
            f1.Pinta_Pixel(-i, y, 1, 1, 1);
        }
        
        
        
    }
    GLvoid drawing_stuff() {

  
            m = (c2[1] - c1[1]) / (c2[0] - c1[0]);
            for (float i = c1[0]; i <= c2[0]; i += 0.05) {
                y = m *( i + -c1[0])+c1[1];
                f1.Pinta_Pixel(i, y, 1, 1, 1);
            }
        
        

    }
    GLvoid DDA() {

        dx = c2[0] - c1[0];
        dy = c2[1] - c1[1];
        if (abs(dx) >= abs(dy))
        steps = abs(dx);
        else
            steps = abs(dy);

        Xinc =(dx / steps);
        Yinc =(dy / steps);
        for (float i = 0; i < steps; i += 0.1) {

            x = c1[0] + Xinc * i;
            y = c1[1] + Yinc * i;
            f1.Pinta_Pixel(x, y, 128, 0, 128);
        }
    }
    GLvoid bresenham() {

        y = c1[1];
        x = c1[0];
        dx = c2[0] - c1[0];
        dy = c2[1] - c1[1];
        float pk = 2 * dy - dx;
        float p1 = 2 * dx - dy;
        int i = 0;
        
        if (abs(dx) >= abs(dy)) {
            steps = abs(dx);
        }
        else {
            steps = abs(dy);
        }
        Xinc = (dx / steps);
        Yinc = (dy / steps);
        if (abs(dx) >= abs(dy)) {
            while (i <= steps) {

                if (pn < 0) {
                    x = x + Xinc;
                    y = y + Yinc;
                    pn = pk + 2 * dy - 2 * dx;
                }
                else {

                    x += Xinc;
                    pn = pk + 2 * dy;

                }
                f1.Pinta_Pixel(x, y, 128, 0, 128);
                i++;
            }
        }
        if (abs(dx) < abs(dy)) {
            while (i <= steps) {

                if (pn > 0) {
                    x = x + Xinc;
                    y = y + Yinc;
                    pn = p1 + 2 * dx - 2 * dy;
                }
                else {

                    y += Yinc;
                    pn = p1 + 2 * dx;

                }
                f1.Pinta_Pixel(x, y, 128, 0, 128);
                i++;
            }
        }

        

    }
    GLvoid setc(int x1, int y1)
    {
        c1[0] = x1;
        c1[1] = y1;
    }
    GLvoid setc2(int x2, int y2)
    {
        c2[0] = x2;
        c2[1] = y2;
    }
};