using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace OptimizationLab1
{
    internal class Program
    {
        #region Main

        private static void Main(string[] args)
        {
            var x0 = new DenseVector(2)
            {
                [0] = -1,
                [1] = 0
            };
            var eps = Math.Pow(10, -5);
            var min = GradientDescentSplitStep(x0, eps);
            
            Console.WriteLine(F(min));
            Console.ReadKey();
        }

        #endregion

        private static double F(double x, double y) //мой вариант
        {
            return x * x + 18 * y * y + 0.01 * x * y + x - y;
        }

        public static double F1(double x, double y)
        {
            return Math.Pow(x, 2) + 8 * Math.Pow(y, 2) + 5 * y - 3;
        }

        private static double F(Vector<double> x)
        {
            return F(x[0], x[1]);
        }


        private static Vector<double> GradientDescentGoldenRatio(Vector<double> x0, double eps)
        {
            var xCur = x0;
            Vector<double> grad = new DenseVector(2)
            {
                [0] = Df_dx(x0[1], x0[1], eps), //using vector x0 fill the grad by-coordinate
                [1] = Df_dy(x0[0], x0[1], eps)
            };

            var alpha = GoldenRatioAlpha(x0, eps, grad / grad.L2Norm()); //eps/10
            var xNext = xCur - grad / grad.L2Norm() * alpha; //нормируем направление шага
            var iteration = 1;


            do
            {
                #region Output

                Console.WriteLine($"Iter {iteration++}: ({xCur[0]}, {xCur[1]}) " + $"Grad: ({grad[0]}, {grad[1]}) " +
                                  $"||Grad||:{grad.L2Norm()} " +
                                  $"||xCur-xNext||:{(xCur - xNext).L2Norm()} " +
                                  $"|F(xCur) - F(xNext)|: {Math.Abs(F(xCur) - F(xNext))}");

                #endregion

                xCur = xNext;

                grad[0] = Df_dx(xCur[0], xCur[1], eps / 10);
                grad[1] = Df_dy(xCur[0], xCur[1], eps / 10);

                alpha = GoldenRatioAlpha(xCur, eps / 10, grad / grad.L2Norm()); //eps/10 + нормируем направление шага

                xNext = xCur - grad / grad.L2Norm() *
                        alpha; //вводим нормированный градиент; ненормированный используем для условия выхода
            } while (!((xCur - xNext).L2Norm() < eps && Math.Abs(F(xCur) - F(xNext)) < eps && grad.L2Norm() < eps));

            Console.WriteLine($"Iter {iteration}: ({xCur[0]}, {xCur[1]}) " + $"Grad: ({grad[0]}, {grad[1]}) " +
                              $"||Grad||:{grad.L2Norm()} " +
                              $"||xCur-xNext||:{(xCur - xNext).L2Norm()} " +
                              $"|F(xCur) - F(xNext)|: {Math.Abs(F(xCur) - F(xNext))}");
            return xNext;
        }
        private static Vector<double> GradientDescentSplitStep(Vector<double> x0, double eps)
        {
            Vector<Double> xCur = x0;
            Vector<double> grad = new DenseVector(2)
            {
                [0] = Df_dx(x0[1], x0[1], eps / 10), //using vector x0 fill the grad by-coordinate
                [1] = Df_dy(x0[0], x0[1], eps / 10)
            };

            double alpha = 0.8*eps; //eps/10
            var xNext = xCur - grad / grad.L2Norm() * alpha; //нормируем направление шага
            int iteration = 1;


            do
            {
                #region Output
                Console.WriteLine($"Iter {iteration++}: ({xCur[0]}, {xCur[1]}) " + $"Grad: ({grad[0]}, {grad[1]}) " +
                                  $"||Grad||:{grad.L2Norm()} " +
                                  $"||xCur-xNext||:{(xCur - xNext).L2Norm()} " +
                                  $"|F(xCur) - F(xNext)|: {Math.Abs(F(xCur) - F(xNext))}");
                #endregion

             
                xCur = xNext;

                grad[0] = Df_dx(xCur[0], xCur[1], eps / 10);
                grad[1] = Df_dy(xCur[0], xCur[1], eps / 10);


                xNext = xCur - grad / grad.L2Norm() *
                        alpha; //вводим нормированный градиент; ненормированный используем для условия выхода

               
                
            } while (!((xCur - xNext).L2Norm() < eps && Math.Abs(F(xCur) - F(xNext)) < eps && grad.L2Norm() < eps));
            Console.WriteLine($"Iter {iteration}: ({xCur[0]}, {xCur[1]}) " + $"Grad: ({grad[0]}, {grad[1]}) " +
                              $"||Grad||:{grad.L2Norm()} " +
                              $"||xCur-xNext||:{(xCur - xNext).L2Norm()} " +
                              $"|F(xCur) - F(xNext)|: {Math.Abs(F(xCur) - F(xNext))}");
            return xNext;
        }

        private static double Df_dx(double x, double y, double h)
        {
            return (F(x + h, y) - F(x - h, y)) / (2 * h);
        }

        private static double Df_dy(double x, double y, double h)
        {
            return (F(x, y + h) - F(x, y - h)) / (2 * h);
        }

        private static double GoldenRatioAlpha(Vector<double> x, double eps, Vector<double> grad)
        {
            double a = 0;
            double b = 1;
            const double G = 1.618d; //constant of a golden ratio
            while (Math.Abs(b - a) > eps)
            {
                var alpha1 = b - (b - a) / G;
                var alpha2 = a + (b - a) / G;
                if (F(x - alpha1 * grad) >= F(x - alpha2 * grad))
                    a = alpha1;
                else
                    b = alpha2;
            }

            return (a + b) / 2;
        }
    }
}