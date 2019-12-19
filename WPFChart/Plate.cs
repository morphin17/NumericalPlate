using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CulculatePlate
{
    public class Plate
    {

        public double a;
        public double b;

        public double h1;
        public double h2;

        public double MaxZ;
        public double MinZ;

        public double tau;
        private double curT;

        public double sigma1;
        public double sigma2;
        private double gamma1;
        private double gamma2;
        public double maxR = 0;


        private Func<double, double> mu;
        private Func<double, double, double> phi;
        private Func<double, double, double, double> f;
        /// <summary>
        /// Пластина в текущий момент времени
        /// </summary>
        public double[,] U;
        public double[,] solutionU;
        public double[,] FF;
        /// <summary>
        /// История развития пластины
        /// </summary>
        public List<double[,]> histiory;

        /// <summary>
        /// Инициализация входных параметров
        /// </summary>
        /// <param name="a">Длина пластины по x</param>
        /// <param name="b">Длина пластины по y</param>
        /// <param name="h1">Шаг по x</param>
        /// <param name="h2">Шаг по y</param>
        /// <param name="tau">Шаг по времени</param>
        /// <param name="phi">Функция при краевом условии t=0</param>
        /// <param name="f">Функция нагрева</param>
        /// <param name="mu">Функция при кравеаом условии на границах пластины</param>
        public Plate(double a, double b, double h1, double h2, double tau,
            Func<double, double, double> phi, Func<double, double, double, double> f, Func<double, double> mu)
        {
            this.a = a;
            this.b = b;
            this.h1 = h1;
            this.h2 = h2;
            this.tau = tau;

            this.sigma1 = 0.5 * (1 - h1 * h1 / (6 * tau));
            this.sigma2 = 0.5 * (1 - h2 * h2 / (6 * tau));
            this.gamma1 = tau / (h1 * h1);
            this.gamma2 = tau / (h2 * h2);

            this.mu = mu;
            this.phi = phi;
            this.f = f;

            this.curT = 0;
            this.U = new double[Convert.ToInt64(a / h1) + 1, Convert.ToInt64(b / h2) + 1];
            this.solutionU = new double[Convert.ToInt64(a / h1) + 1, Convert.ToInt64(b / h2) + 1];
            this.FF = new double[Convert.ToInt64(a / h1) + 1, Convert.ToInt64(b / h2) + 1];

            this.MaxZ = 0;
            this.MinZ = Double.MaxValue;
            //Пластина в нулевой момент времени
            for (int i = 0; i < this.U.GetLength(0); i++)
            {
                for (int j = 0; j < this.U.GetLength(1); j++)
                {
                    this.U[i, j] = phi(h1 * i, h2 * j);
                }
            }
            
            if (this.MaxZ == 0)
            {
                this.MaxZ = 1;
            }

            //Краевые условия на краях пластины
            this.conditionT();

            for (int i = 0; i < this.U.GetLength(0); i++)
            {
                for (int j = 0; j < this.U.GetLength(1); j++)
                {
                    if (this.MaxZ < this.U[i, j])
                    {
                        this.MaxZ = this.U[i, j];
                    }

                    if (this.MinZ > this.U[i, j])
                    {
                        this.MinZ = this.U[i, j];
                    }
                }
            }

            this.histiory = new List<double[,]>();
            this.histiory.Add(Copy(this.U));
        }

        private void calcFF()
        {
            for (int i = 0; i < this.U.GetLength(0); i++)
            {
                for (int j = 0; j < this.U.GetLength(1); j++)
                {
                    this.FF[i, j] = f(h1 * i, h2 * j, curT + tau/2);
                }
            }
        }
        /// <summary>
        /// Посчитать сетку для аналитического решения
        /// </summary>
        /// <param name="solution"></param>
        public void calcSolution(Func<double, double, double, double> solution)
        {
            for (int i = 0; i < this.U.GetLength(0); i++)
            {
                for (int j = 0; j < this.U.GetLength(1); j++)
                {
                    this.solutionU[i, j] = solution(h1 * i, h2 * j, curT);
                }
            }
        } 

        /// <summary>
        /// Краевые условия на краях пластины
        /// </summary>
        private void conditionT()
        {
            for (int i = 0; i < U.GetLength(0); i++)
            {
                U[i, 0] = mu(curT);
                U[i, U.GetLength(1) - 1] = mu(curT);
            }

            for (int j = 0; j < U.GetLength(1); j++)
            {
                U[0, j] = mu(curT);
                U[U.GetLength(0) - 1, j] = mu(curT);
            }
        }

        /// <summary>
        /// Вспомогательная функция для подсчета схемы
        /// </summary>
        /// <param name="U">Пластина</param>
        /// <param name="i">Координаты</param>
        /// <param name="j">Координаты</param>
        /// <returns></returns>
        private double F(double[,] U, int i, int j)
        {
            //return (1 - sigma1) * (U[i - 1, j] - 2 * U[i, j] + U[i + 1, j]) / (h1 * h1) +
            //    (1 - sigma2) * (U[i, j - 1] - 2 * U[i, j] + U[i, j + 1]) / (h2 * h2) +
            //    (1 - sigma1) * (1 - sigma2) * tau *
            //    (U[i - 1, j - 1] - 2 * U[i - 1, j] + U[i - 1, j + 1] - 2 * U[i, j - 1] + 4 * U[i, j] - 2 * U[i, j + 1] + U[i + 1, j - 1] - 2 * U[i + 1, j] + U[i + 1, j + 1]) /
            //    ((h1 * h1) * (h2 * h2)) +
            //    ((double)1 / 12) * (FF[i - 1, j] - 4 * FF[i, j] + FF[i + 1, j] + FF[i, j - 1] + FF[i, j + 1]) + FF[i, j];
            return (1 - sigma1) * (U[i - 1, j] - 2 * U[i, j] + U[i + 1, j]) / (h1 * h1) +
                (1 - sigma2) * (U[i, j - 1] - 2 * U[i, j] + U[i, j + 1]) / (h2 * h2) +
                (1 - sigma1)*(1 - sigma2) * tau *
                (U[i - 1, j - 1] - 2 * U[i - 1, j] + U[i - 1, j + 1] - 2 * U[i, j - 1] + 4 * U[i, j] - 2 * U[i, j + 1] + U[i + 1, j - 1] - 2 * U[i + 1, j] + U[i + 1, j + 1]) /
                ((h1 * h1) * (h2 * h2)) +
                ((double)1 / 12) * (FF[i - 1, j] - 4 * FF[i, j] + FF[i + 1, j] + FF[i, j - 1] + FF[i, j + 1]) + FF[i, j];
        }

        /// <summary>
        /// Переход на следующий временной слой
        /// </summary>
        public void nextU()
        {
            calcFF();
            Double[,] N;
            double[] alpha;
            double[] beta;

            N = Copy(U);
            alpha = new double[U.GetLength(0) + 1];
            beta = new double[U.GetLength(0) + 1];

            //Прогонка по x
            for (int j = 1; j < U.GetLength(1) - 1; j++)
            {
                alpha[1] = 0;
                beta[1] = U[0, j];

                for (int i = 1; i < U.GetLength(0) - 1; i++)
                {
                    alpha[i + 1] = sigma1 * gamma1 / (1 + 2 * sigma1 * gamma1 - alpha[i] * sigma1 * gamma1);
                    beta[i + 1] = (tau * F(N, i, j) + N[i, j] + sigma1 * gamma1 * beta[i]) / (1 + 2 * sigma1 * gamma1 - alpha[i] * sigma1 * gamma1);
                }

                for (int i = U.GetLength(0) - 2; i > 0; i--)
                {
                    U[i, j] = alpha[i + 1] * U[i + 1, j] + beta[i + 1];
                }
            }

            N = Copy(U);
            alpha = new double[U.GetLength(1) + 1];
            beta = new double[U.GetLength(1) + 1];

            //Прогонка по y
            for (int i = 1; i < U.GetLength(0) - 1; i++)
            {
                alpha[1] = 0;
                beta[1] = U[i, 0];

                for (int j = 1; j < U.GetLength(1) - 1; j++)
                {
                    alpha[j + 1] = sigma2 * gamma2 / (1 + 2 * sigma2 * gamma2 - alpha[j] * sigma2 * gamma2);
                    beta[j + 1] = ( N[i, j] + sigma2 * gamma2 * beta[j]) / (1 + 2 * sigma2 * gamma2 - alpha[j] * sigma2 * gamma2);
                }

                for (int j = U.GetLength(1) - 2; j > 0; j--)
                {
                    U[i, j] = alpha[j + 1] * U[i, j + 1] + beta[j + 1];
                }
            }

            bool flag = false;
            for (int i = 0; i < this.U.GetLength(0); i++)
            {
                for (int j = 0; j < this.U.GetLength(1); j++)
                {
                    if (this.MaxZ < this.U[i, j])
                    {
                        flag = true;
                        break;
                    }
                }
            }

            if (flag)
            {
                this.MaxZ += 3;
            }

            curT += tau;
            conditionT();
            histiory.Add(Copy(U));
        }

        
        public double errorValue()
        {
            double error = 0;

            for (int i = 0; i < U.GetLength(0); i++)
            {
                for (int j = 0; j < U.GetLength(1); j++)
                {
                    if (Math.Abs(U[i, j] - solutionU[i,j]) > error)
                    {
                        error = Math.Abs(U[i, j] - solutionU[i, j]);
                    }
                }
            }

            if (maxR < error)
            {
                maxR = error;
            }

            return error;
        }

        static T[,] Copy<T>(T[,] array)
        {
            T[,] newArray = new T[array.GetLength(0), array.GetLength(1)];
            for (int i = 0; i < array.GetLength(0); i++)
                for (int j = 0; j < array.GetLength(1); j++)
                    newArray[i, j] = array[i, j];
            return newArray;
        }

    }

}