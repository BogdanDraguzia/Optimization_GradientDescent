using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OptimizationLab1;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using  lab1 = OptimizationLab1.Program;
namespace OptimizationLab3_GradientProjection
{
    class Program
    {
        private static double F(double x, double y, double z) => 4 * Math.Pow(x, 2) + Math.Pow(y, 2) + Math.Pow(z, 2);

        private static double F(Vector<double> xVector)
        {
            return F(xVector[0], xVector[1], xVector[2]);
        }
       
        
        //(p,a)=beta
        public static readonly Vector<double> HyperPlane = new DenseVector(new Double[] {1, 2, 3});
        

        
        public static Vector<double> GradientProjection(Vector<double> x0, double eps)
        {
            var xCur = x0.GetProjectionHyperplane(HyperPlane,1);
            Vector<double> grad = new DenseVector(3)
            {
                [0] = Df_dx(x0[0], x0[1], x0[2], eps), //using vector x0 fill the grad by-coordinate
                [1] = Df_dy(x0[0], x0[1], x0[2], eps),
                [2] = Df_dz(x0[0], x0[1], x0[2], eps)
            };


            //подбираем длину шага в направлении проекции
            ////var alpha = GoldenRatioAlpha(x0, eps, grad.Normalize(2).GetProjectionHyperplane(HyperPlane,1)); //project???
            //.GetProjectionHyperplane(HyperPlane, 1)

            var alpha = 1 / 2;

            //(HyperPlane,a)=1
            var xNext = (xCur - alpha*grad).GetProjectionHyperplane(HyperPlane, 1); //нормируем направление шага
            var iteration = 1;//not normalized


            do
            {
                #region Output

                Console.WriteLine($"Iter {iteration++}: xCur:{xCur.ToStR()}" + $"F = {F(xCur)}" +
                                  $"Grad:{(-grad).ToStR()} " +
                                  $"||Grad||:{grad.L2Norm()} " +
                                  $"||xCur-xNext||:{(xCur - xNext).L2Norm()}");
                #endregion

                if (!IsOnPlane(HyperPlane, xCur, eps))
                {
                    throw new Exception();
                    
                }
                xCur = xNext;

                grad[0] = Df_dx(xCur[0], xCur[1], xCur[2], eps / 10);
                grad[1] = Df_dy(xCur[0], xCur[1], xCur[2], eps / 10);
                grad[2] = Df_dz(xCur[0], xCur[1], xCur[2], eps / 10);

                //подбираем длину шага в направлении проекции
                ////alpha = GoldenRatioAlpha(xCur, eps / 10, grad.Normalize(2).GetProjectionHyperplane(HyperPlane, 1)); //eps/10 + нормируем направление шага
                //.GetProjectionHyperplane(HyperPlane,1)

                alpha = 1 /(iteration+1);

                xNext = (xCur - alpha*grad).GetProjectionHyperplane(HyperPlane, 1); //вводим нормированный градиент; ненормированный используем для условия выхода
            }while (!((xCur - xNext).L2Norm() < eps)); //not normalized now
            
            Console.WriteLine($"Iter {iteration}: xCur:{xCur.ToStR()} "+
                              $"||xCur-xNext||:{(xCur - xNext).L2Norm()} " +
                              $"|F(xCur) - F(xNext)|: {Math.Abs(F(xCur) - F(xNext))}");
            return xNext;
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
        
        private static double Df_dx(double x, double y, double z, double h)
        {
            return (F(x + h, y, z) - F(x - h, y, z)) / (2 * h);
        }
        private static double Df_dy(double x, double y, double z, double h)
        {
            return (F(x, y + h, z) - F(x, y - h, z)) / (2 * h);
        }
        private static double Df_dz(double x, double y, double z, double h)
        {
            return (F(x, y, z + h) - F(x, y, z - h)) / (2 * h);
        }

        public static bool IsOnPlane(Vector<double> plane,Vector<double> x, double eps)
        {
            return Math.Abs(plane * x - 1) < eps;
        }


        static void Main(string[] args)
        {
            Vector<double> x0 = new DenseVector(3)
            {
                [0] = 2,
                [1] = 2, 
                [2] = 2
            };
            var eps = Math.Pow(10, -5);
            GradientProjection(x0, eps);
           

            Console.ReadKey();
        }
    }

    public static class VectorString //extension method for vector output
    {
        public static string ToStR(this Vector<double> x) => $"({x[0]}, {x[1]}, {x[2]})";

        public static Vector<double> GetProjectionHyperplane(this Vector<double> a,Vector<double> p, double beta)
        {
            return a + (beta - p * a) * p / Math.Pow(p.L2Norm(), 2);
        }
    }

}

