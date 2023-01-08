#include "problem/wfg.h"

#include <cmath>
#include <cstdlib>

#include "core/global.h"

namespace emoc
{
    /**
     * transformation functions.
     * Shift: Linear
     * @param y
     * @param A
     * @return
     */
    double s_linear(double y, double A)
    {
        double result = fabs(y - A) / fabs(floor(A - y) + A);
        if (result > 1) result = 1;
        if (result < 0) result = 0;
        return result;
    }

    /**
     *
     * @param y primary parameters
     * @param y_size size of primary parameters (k+l)
     * @param k number of position-related parameters
     * @param result t^1
     * @return y_size size of primary parameters (k+l)
     */
    int WFG1_t1(double *y, int y_size, int k, double *result)
    {
        int i;
        for (i = 0; i < k; i++)
            result[i] = y[i];

        for (i = k; i < y_size; i++)
            result[i] = s_linear(y[i], 0.35);

        return y_size;
    }
    /**
     * y -> t1
     * @param y origin y
     * @param y_length k+l
     * @param k number of position-related parameters
     * @param t1 save the results
     */
    void y_2_t1(double *y, const int &y_length, const int &k, double *t1)
    {
        int i=0;
        for (; i < k; ++i)
        {
            t1[i] = y[i];
        }
        for (; i < y_length; ++i)
        {
            t1[i] = s_linear(y[i], 0.35);
        }
    }
    /**
     * transformation function.
     * Reduction: Non-separable
     * @param y y_{k+2(i-k)-2}
     * @param size 2
     * @param A 2
     * @return
     */
    double r_nonsep(double *y, int size, const int A)
    {
        int i, j;
        double result;
        double numerator = 0.0;
        for (i = 0; i < size; i++)
        {
            numerator += y[i];

            for (j = 0; j <= A - 2; j++)
                numerator += fabs(y[i] - y[(i + j + 1) % size]);
        }

        const double tmp = ceil(A / 2.0);
        const double denominator = size * tmp * (1.0 + 2.0 * A - 2.0 * tmp) / A;

        result = numerator / denominator;
        if (result > 1) result = 1;
        if (result < 0) result = 0;

        return result;
    }

    /**
     *
     * @param y t^1
     * @param y_size n:=k+l (dec_num_)
     * @param k number of position-related parameters
     * @param result t^2
     * @return (k + l / 2)
     */
    int WFG2_t2(double *y, int y_size, int k, double *result)
    {
        int i;
        const int l = y_size - k;// number of distance-related parameters

        for (i = 0; i < k; i++)
            result[i] = y[i];

        for (i = k + 1; i <= k + l / 2; i++)
        {
            const int head = k + 2 * (i - k) - 2;
            result[i] = r_nonsep(y + head, 2, 2);//TODO: DEBUG result[i-1]
        }

        return k + l / 2;//TODO: DEBUG
    }
    void t1_2_t2(double *t1, const int &t1_length, const int &k, const int &l, double *t2)
    {
        int i = 0;
        for (; i < k; ++i)
        {
            t2[i] = t1[i];
        }
        for (; i < k+l/2; ++i)
        {
            const int head = k + 2 * (i + 1 - k) - 2;

            t2[i] = r_nonsep(t1 + head, 2, 2);
        }
    }

    /**
     *
     * @param y y_k
     * @param y_size l/2
     * @param w
     * @return
     */
    double r_sum(double *y, int y_size)
    {
        int i;
        double result;
        double numerator = 0.0;
        double denominator = 0.0;

        for (i = 0; i < y_size; i++)
        {
            numerator += y[i];
            denominator += 1.0;
        }

        result = numerator / denominator;
        if (result > 1) result = 1;
        if (result < 0) result = 0;

        return result;
    }

    /**
     *
     * @param y t^2
     * @param y_size k + l / 2
     * @param k
     * @param M nobj
     * @param result t^3
     * @param wfg_w
     * @return M:=nobj
     */
    int WFG2_t3(double *y, int y_size, int k, int M, double *result, double *wfg_w)
    {
        int i;

        for (i = 0; i < y_size; i++)
            wfg_w[i] = 1.0;

        for (i = 1; i <= M - 1; i++)
        {
            const int head = (i - 1) * k / (M - 1);
            const int tail = i * k / (M - 1);
            result[i - 1] = r_sum(y + head, tail - head);
        }

        result[M - 1] = r_sum(y + k, y_size - k);

        return M;
    }
    void t2_2_t3(double *t2, const int &k, const int &l, const int &M, double *t3)
    {
        for (int i = 0; i < M-1; ++i)
        {
            const int head = i * k / (M - 1);
            const int tail = (i+1) * k / (M - 1);
            t3[i] = r_sum(t2 + head, tail - head);
        }
        t3[M-1] = r_sum(t2+k,l/2);
    }

    /**
     *
     * @param x t^3 -> transition vector
     * @param result x_{1,...,M} -> underlying vector
     * @param size nobj -> the number of objectives M
     * @param Degenerate 1 -> whether to degenerate
     */
    void calculate_x(double *x, double *result, int size, int Degenerate)
    {
        int i;
        double val = x[size - 1];

        if (!Degenerate && val < 1)
        {
            val = 1;
        }

        result[0] = x[0];
        result[size - 1] = x[size - 1];

        for (i = 1; i < size - 1; i++)
        {
            result[i] = (x[i] - 0.5) * val + 0.5;
        }
    }

    /**
     * Shape function.
     * Linear
     * @param x x_{1,...,m-1}
     * @param M nobj
     * @param m
     * @return h_i
     */
    double linear(double *x, int M, int m)
    {
        int i;
        double result = 1.0;

        for (i = 1; i <= M - m; i++)
            result *= x[i - 1];

        if (m != 1)
            result *= 1 - x[M - m];

        // FIXME
        if (result > 1) result = 1;
        if (result < 0) result = 0;

        return result;
    }

    /**
     *
     * @param D 1.0
     * @param x x_M
     * @param h obj
     * @param size nobj
     * @param result obj
     */
    void calculate_f(double D, double x, double *h, int size, double *result)
    {
        int i;
        int S = 0;

        for (i = 0; i < size; i++)
        {
            S = S + 2;
            result[i] = D * x + S * h[i];
        }
    }

    /**
     *
     * @param y
     * @param y_size nobj
     * @param result obj
     * @param temp
     * @param Degenerate
     */
    void WFG3_shape(double *y, int y_size, double *result, double *temp, int Degenerate)
    {
        int i;

        calculate_x(y, temp, y_size, Degenerate);

        for (i = 1; i <= y_size; i++)
            result[i - 1] = linear(temp, y_size, i);

        calculate_f(1.0, temp[y_size - 1], result, y_size, result);
    }
    void WFG3_final(double *t3, const int &M, int Degenerate, double *final_obj)
    {
        double x[M];
        calculate_x(t3, x, M, Degenerate);
        for (int i = 0; i < M; ++i)
        {
            //PF :sum of current final_obj is 1.
            final_obj[i] = linear(x, M, i+1);
        }
        calculate_f(1.0, x[M-1], final_obj, M, final_obj);
    }

    WFG3::WFG3(int dec_num, int obj_num) : Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    WFG3::~WFG3()
    {
    }
    /**
     * calculate objective value of an individual for the wfg3.
     * n:=dec_num_;
     * k:=2(M-1);
     * l:=20;
     * k+l=n;
     * @param ind
     */
    void WFG3::CalObj(Individual *ind)
    {
        // A_1=1, A_{2,...,M-1}=0
        const int Degenerate = 1;
        // k: length of position-related parameters.
        const int k = 2*(obj_num_-1);
        // l: length of distance-related parameters, even required.
        const int l = dec_num_ - k;
        const int t1_length = dec_num_;
        const int t2_length = k + l / 2;
        const int t3_length = obj_num_;
        double t1[t1_length];
        double t2[t2_length];
        double t3[t3_length];

        // n = k + l -> 2 * (M - 1) + 20
        double *xreal = ind->dec_;
        const int y_length = dec_num_;
        double *obj = ind->obj_;

        y_2_t1(xreal, y_length, k, t1);
        t1_2_t2(t1, t1_length, k, l, t2);
        t2_2_t3(t2, k, l, t3_length, t3);
        WFG3_final(t3, obj_num_, Degenerate, obj);
    }

    WFG1::WFG1(int dec_num, int obj_num) : Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    WFG1::~WFG1()
    {
    }

    void WFG1::CalObj(Individual *ind)
    {
    }

    WFG2::WFG2(int dec_num, int obj_num) : Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    WFG2::~WFG2()
    {
    }

    void WFG2::CalObj(Individual *ind)
    {
    }

    WFG4::WFG4(int dec_num, int obj_num) : Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    WFG4::~WFG4()
    {
    }

    void WFG4::CalObj(Individual *ind)
    {
    }

    WFG5::WFG5(int dec_num, int obj_num) : Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    WFG5::~WFG5()
    {
    }

    void WFG5::CalObj(Individual *ind)
    {
    }

    WFG6::WFG6(int dec_num, int obj_num) : Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    WFG6::~WFG6()
    {
    }

    void WFG6::CalObj(Individual *ind)
    {
    }

    WFG7::WFG7(int dec_num, int obj_num) : Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    WFG7::~WFG7()
    {
    }

    void WFG7::CalObj(Individual *ind)
    {
    }

    WFG8::WFG8(int dec_num, int obj_num) : Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    WFG8::~WFG8()
    {
    }

    void WFG8::CalObj(Individual *ind)
    {
    }

    WFG9::WFG9(int dec_num, int obj_num) : Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    WFG9::~WFG9()
    {
    }

    void WFG9::CalObj(Individual *ind)
    {
    }
}