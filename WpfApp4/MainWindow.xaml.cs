using System;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Media3D;
using CulculatePlate;
using System.Windows.Threading;
using System.Collections;
using WPFChart3D;
using LiveCharts;
using LiveCharts.Wpf;

namespace WpfApp4
{

    public class Chart3D
    {
        // transform class object for rotate the 3d model
        public WPFChart3D.TransformMatrix m_transformMatrix = new WPFChart3D.TransformMatrix();

        // ***************************** 3d chart ***************************
        private WPFChart3D.Chart3D m_3dChart;       // data for 3d chart
        public int m_nChartModelIndex = -1;         // model index in the Viewport3d
        Viewport3D mainViewport;

        // ***************************** selection rect ***************************
        ViewportRect m_selectRect = new ViewportRect();
        public int m_nRectModelIndex = -1;

        public Chart3D(Viewport3D mainViewport)
        {
            this.mainViewport = mainViewport;
            // selection rect
            m_selectRect.SetRect(new Point(-0.5, -0.5), new Point(-0.5, -0.5));
            WPFChart3D.Model3D model3d = new WPFChart3D.Model3D();
            ArrayList meshs = m_selectRect.GetMeshes();
            m_nRectModelIndex = model3d.UpdateModel(meshs, null, m_nRectModelIndex, this.mainViewport);
        }

        public void plotNumerical(Plate p)
        {
            int nXNo = p.U.GetLength(0);
            int nYNo = p.U.GetLength(1);

            //Устанавливаем сетку
            m_3dChart = new UniformSurfaceChart3D();
            ((UniformSurfaceChart3D)m_3dChart).SetGrid(nXNo, nYNo, 0, (float)p.a, 0, (float)p.b, p.U);
            m_3dChart.GetDataRange();

            int nVertNo = m_3dChart.GetDataNo();
            double zMin = p.MinZ;
            double zMax = p.MaxZ;
            for (int i = 0; i < nVertNo; i++)
            {
                Vertex3D vert = m_3dChart[i];
                double h = (vert.z - zMin) / (zMax - zMin);

                Color color = WPFChart3D.TextureMapping.PseudoColor(h);
                m_3dChart[i].color = color;
            }


            ArrayList meshs = ((UniformSurfaceChart3D)m_3dChart).GetMeshes();

            // 5. display vertex no and triangle no of this surface chart
            UpdateModelSizeInfo(meshs);

            // 6. Set the model display of surface chart
            WPFChart3D.Model3D model3d = new WPFChart3D.Model3D();
            Material backMaterial = new DiffuseMaterial(new SolidColorBrush(Colors.Gray));
            m_nChartModelIndex = model3d.UpdateModel(meshs, backMaterial, m_nChartModelIndex, this.mainViewport);

            // 7. set projection matrix, so the data is in the display region
            float xMin = m_3dChart.XMin();
            float xMax = m_3dChart.XMax();
            m_transformMatrix.CalculateProjectionMatrix(xMin, xMax, xMin, xMax, zMin, zMax, 0.5);
            TransformChart();
        }

        public void plotAnalytic(Plate p)
        {
            int nXNo = p.U.GetLength(0);
            int nYNo = p.U.GetLength(1);

            //Устанавливаем сетку
            m_3dChart = new UniformSurfaceChart3D();
            ((UniformSurfaceChart3D)m_3dChart).SetGrid(nXNo, nYNo, 0, (float)p.a, 0, (float)p.b, p.solutionU);
            m_3dChart.GetDataRange();

            int nVertNo = m_3dChart.GetDataNo();
            double zMin = p.MinZ;
            double zMax = p.MaxZ;
            for (int i = 0; i < nVertNo; i++)
            {
                Vertex3D vert = m_3dChart[i];
                double h = (vert.z - zMin) / (zMax - zMin);

                Color color = WPFChart3D.TextureMapping.PseudoColor(h);
                m_3dChart[i].color = color;
            }


            ArrayList meshs = ((UniformSurfaceChart3D)m_3dChart).GetMeshes();

            // 5. display vertex no and triangle no of this surface chart
            UpdateModelSizeInfo(meshs);

            // 6. Set the model display of surface chart
            WPFChart3D.Model3D model3d = new WPFChart3D.Model3D();
            Material backMaterial = new DiffuseMaterial(new SolidColorBrush(Colors.Gray));
            m_nChartModelIndex = model3d.UpdateModel(meshs, backMaterial, m_nChartModelIndex, this.mainViewport);

            // 7. set projection matrix, so the data is in the display region
            float xMin = m_3dChart.XMin();
            float xMax = m_3dChart.XMax();
            m_transformMatrix.CalculateProjectionMatrix(xMin, xMax, xMin, xMax, zMin, zMax, 0.5);
            TransformChart();
        }


        private void UpdateModelSizeInfo(ArrayList meshs)
        {
            int nMeshNo = meshs.Count;
            int nChartVertNo = 0;
            int nChartTriangelNo = 0;
            for (int i = 0; i < nMeshNo; i++)
            {
                nChartVertNo += ((Mesh3D)meshs[i]).GetVertexNo();
                nChartTriangelNo += ((Mesh3D)meshs[i]).GetTriangleNo();
            }
            //labelVertNo.Content = String.Format("Vertex No: {0:d}", nChartVertNo);
            //labelTriNo.Content = String.Format("Triangle No: {0:d}", nChartTriangelNo);
        }

        // this function is used to rotate, drag and zoom the 3d chart
        public void TransformChart()
        {
            if (m_nChartModelIndex == -1) return;
            ModelVisual3D visual3d = (ModelVisual3D)(this.mainViewport.Children[m_nChartModelIndex]);
            if (visual3d.Content == null) return;
            Transform3DGroup group1 = visual3d.Content.Transform as Transform3DGroup;
            group1.Children.Clear();
            group1.Children.Add(new MatrixTransform3D(m_transformMatrix.m_totalMatrix));
        }
    }

    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        

        /** Анимация **/
        private bool animateFlag = false;       //Воспроизвестиъ
        private bool solutionFlag = false;      //Дано ли аналитическое решение
        DispatcherTimer timer;                  //Таймер для анимации
        double t = 0;                           //Время
        double runge = 0;


        public Plate p;                         //Пластина, числнное решение
        public Plate p2;                        //Пластина, для шага в два раза больше
        

        Chart3D chart1;                         //График - Аналитическое решение
        Chart3D chart2;                         //График - Численное решение
        public SeriesCollection SeriesCollection { get; set; }  //График - Погрешность

        //Спрятать график Анлитического решения и график погрешности
        private void hideChart()
        {
            mainViewport.Visibility = Visibility.Hidden;
            textChart1.Visibility = Visibility.Visible;
            //Graph1.Visibility = Visibility.Hidden;
            //textGraph1.Visibility = Visibility.Visible;
        }

        //Показать 
        private void showChart()
        {
            mainViewport.Visibility = Visibility.Visible;
            textChart1.Visibility = Visibility.Hidden;
            Graph1.Visibility = Visibility.Visible;
            textGraph1.Visibility = Visibility.Hidden;
        }


        public MainWindow()
        {
            InitializeComponent();
            hideChart();

            chart1 = new Chart3D(this.mainViewport);
            chart2 = new Chart3D(this.mainViewport2);

            timer = new DispatcherTimer();

            timer.Tick += new EventHandler(timer_Tick);

            timer.Interval = new TimeSpan(0, 0, 0, 0, 40);

            timer.Start();

            SeriesCollection = new SeriesCollection();

            ////YFormatter = value => value.ToString("C");

            ////modifying the series collection will animate and update the chart
            SeriesCollection.Add(new LineSeries
            {
                Title = "Погрешность",
                Values = new ChartValues<double>(),
                PointForeground = Brushes.Gray,
            });

            SeriesCollection.Add(new LineSeries
            {
                Title = "Погрешность по Рунге",
                Values = new ChartValues<double>(),
                PointForeground = Brushes.Gray
            });
            ////modifying any series values will also animate and update the chart
            DataContext = this;
        }

        private void timer_Tick(object sender, EventArgs e)
        {
            
            if (animateFlag)
            {
                try
                {
                    p.nextU();  //Переходим на следущий временной слой по h

                }
                catch
                {
                    MessageBox.Show("Ошибка. Проверьте корректны ли введены начальные условия.");
                    animateFlag = false;
                    if (animateFlag)
                    {
                        Play.Content = "⏸";
                    }
                    else
                    {
                        Play.Content = "▶";
                    }
                    return;
                }
                p2.nextU(); //Переходим на следущий временной слой по 2h

                //Если дано решение, рисуем
                if (solutionFlag)
                {
                    p.calcSolution(solution);     //Считаем решение на сетке

                    chart1.plotAnalytic(p);       //Рисуем
                    chart1.TransformChart();

                    double r = p.errorValue();    //Погрешность численного решение отностительно аналитического  
                    SeriesCollection[0].Values.Add(r); //Добавляем точку на график
                }
                
                //Рисуем численное решение
                chart2.plotNumerical(p);
                chart2.TransformChart();

                //Считаем погрешность по рунге
                double tempR = rungeR(p, p2);
                runge = Math.Max(tempR, runge);
                SeriesCollection[1].Values.Add(tempR);

                //Обнавляем погрешность на текущий момент времени
                maxRText.Text = "R = " + p.maxR.ToString("N6");
                maxRRText.Text = "R2 = " + runge.ToString("N6");
            }

            //Вращение графиков
            //Если дано решение, рисуем
            if (solutionFlag)
            {
                chart1.m_transformMatrix.rotate();
                chart1.TransformChart();
            }

            chart2.m_transformMatrix.rotate();
            chart2.TransformChart();

        }
        
        

        NCalc.Expression exprF;
        private double f(double x, double y, double t)
        {

            exprF = new NCalc.Expression(fText.Text);

            exprF.Parameters["pi"] = Math.PI;
            exprF.Parameters["e"] = Math.E;
            exprF.Parameters["a"] = Convert.ToDouble(aText.Text);
            exprF.Parameters["b"] = Convert.ToDouble(bText.Text);
            exprF.Parameters["x"] = x;
            exprF.Parameters["y"] = y;
            exprF.Parameters["t"] = t;

            double r = Convert.ToDouble(exprF.Evaluate());
            return r;
        }

        NCalc.Expression exprPHI;
        private double phi(double x, double y)
        {
            exprPHI = new NCalc.Expression(phiText.Text);

            exprPHI.Parameters["pi"] = Math.PI;
            exprPHI.Parameters["e"] = Math.E;
            exprPHI.Parameters["a"] = Convert.ToDouble(aText.Text);
            exprPHI.Parameters["b"] = Convert.ToDouble(bText.Text);
            exprPHI.Parameters["x"] = x;
            exprPHI.Parameters["y"] = y;

            double r = Convert.ToDouble(exprPHI.Evaluate());
            return r;
        }

        NCalc.Expression exprMU;
        private double mu(double t)
        {
            exprMU = new NCalc.Expression(muText.Text);
            exprMU.Parameters["pi"] = Math.PI;
            exprMU.Parameters["e"] = Math.E;
            exprMU.Parameters["a"] = Convert.ToDouble(aText.Text);
            exprMU.Parameters["b"] = Convert.ToDouble(bText.Text);
            exprMU.Parameters["t"] = t;

            double r = Convert.ToDouble(exprMU.Evaluate());
            return r;
        }

        NCalc.Expression exprU;
        private double solution(double x, double y, double t)
        {
            
            exprU = new NCalc.Expression(solutionText.Text);
            exprU.Parameters["pi"] = Math.PI;
            exprU.Parameters["e"] = Math.E;
            exprU.Parameters["a"] = Convert.ToDouble(aText.Text);
            exprU.Parameters["b"] = Convert.ToDouble(bText.Text);
            exprU.Parameters["t"] = t;
            exprU.Parameters["x"] = x;
            exprU.Parameters["y"] = y;

            double r = Convert.ToDouble(exprU.Evaluate());
            return r;
        }

        private void build_Click(object sender, RoutedEventArgs e)
        {
            if (solutionText.Text == "")
            {
                solutionFlag = false;
            }
            else
            {
                solutionFlag = true;
            }

            try
            { 
                p = new Plate(
                    Convert.ToDouble(aText.Text),
                    Convert.ToDouble(bText.Text), 
                    Convert.ToDouble(h1Text.Text), 
                    Convert.ToDouble(h2Text.Text), 
                    Convert.ToDouble(tauText.Text), 
                    phi, 
                    f, 
                    mu
                );

                p2 = new Plate(
                    Convert.ToDouble(aText.Text),
                    Convert.ToDouble(bText.Text),
                    Convert.ToDouble(h1Text.Text)*2,
                    Convert.ToDouble(h2Text.Text)*2,
                    Convert.ToDouble(tauText.Text),
                    phi,
                    f,
                    mu
                );

            }
            catch
            {
                MessageBox.Show("Ошибка. Проверьте корректно ли вы ввели начальные условия.");
                return;
            }
            SeriesCollection[0].Values.Clear();
            SeriesCollection[1].Values.Clear();


            //Если дано решение, рисуем
            if (solutionFlag)
            {
                showChart();
                p.calcSolution(solution);     //Считаем решение на сетке

                chart1.plotAnalytic(p);       //Рисуем
                chart1.TransformChart();

                double r = p.errorValue();    //Погрешность численного решение отностительно аналитического  
                SeriesCollection[0].Values.Add(r); //Добавляем точку на график
            }
            else
            {
                hideChart();
            }

            //Рисуем численное решение
            chart2.plotNumerical(p);
            chart2.TransformChart();

            //Считаем погрешность по рунге
            double tempR = rungeR(p, p2);
            runge = Math.Max(tempR, runge);
            SeriesCollection[1].Values.Add(tempR);

            //Обнавляем погрешность на текущий момент времени
            maxRText.Text = "R = " + p.maxR.ToString("N6");
            maxRRText.Text = "R2 = " + runge.ToString("N6");
            sigma1Text.Text = "σ1 = " + p.sigma1.ToString("N6");
            sigma2Text.Text = "σ2 = " + p.sigma2.ToString("N6");
        }

        private double rungeR(Plate p, Plate p2)
        {
            double runge = 0;
            for (int i = 0; i < p.U.GetLength(0) / 2; i += 2)
            {
                for (int j = 0; j < p.U.GetLength(1) / 2; j += 2)
                {
                    double temp = Math.Abs(p2.U[i, j] - p.U[2 * i, 2 * j]);
                    if (temp > runge)
                    {
                        runge = temp;
                    }
                }
            }

            runge = runge*100000* p.tau*p.tau/15;
            return runge;
        }

        private void MenuItem_Test1(object sender, RoutedEventArgs e)
        {
            fText.Text = "0";
            phiText.Text = "Sin(pi*x/a)*Sin(pi*y/b)";
            muText.Text = "0";
            aText.Text = "5";
            bText.Text = "5";
            solutionText.Text = "Sin(pi*x/a)*Sin(pi*y/b)*Pow(e, -(pi*pi/(a*a) + pi*pi/(b*b))*t)";
            h1Text.Text = "0,1";
            h2Text.Text = "0,1";
            tauText.Text = "0,1";
        }

        private void MenuItem_Test2(object sender, RoutedEventArgs e)
        {
            fText.Text = "Sin(pi*x/a)*Sin(pi*y/b)";
            phiText.Text = "0";
            muText.Text = "0";
            aText.Text = "5";
            bText.Text = "5";
            solutionText.Text = "((1/(pi*pi/(a*a) + pi*pi/(b*b)))*Pow(e, (pi*pi/(a*a) + pi*pi/(b*b))*t) - (1/(pi*pi/(a*a) + pi*pi/(b*b)))*Pow(e, (pi*pi/(a*a) + pi*pi/(b*b))*0))*Sin(pi*x/a)*Sin(pi*y/b)*Pow(e, -(pi*pi/(a*a) + pi*pi/(b*b))*t)";
            h1Text.Text = "0,1";
            h2Text.Text = "0,1";
            tauText.Text = "0,1";
        }

        private void MenuItem_Test3(object sender, RoutedEventArgs e)
        {
            fText.Text = "Sin(pi*x/a)*Sin(pi*y/b)*t";
            phiText.Text = "0";
            muText.Text = "0";
            aText.Text = "5";
            bText.Text = "5";
            solutionText.Text = "(Pow(e,(pi*pi/(a*a) + pi*pi/(b*b))*t)*((pi*pi/(a*a) + pi*pi/(b*b))*t -1)/Pow((pi*pi/(a*a) + pi*pi/(b*b)),2)-Pow(e,(pi*pi/(a*a) + pi*pi/(b*b))*0)*((pi*pi/(a*a) + pi*pi/(b*b))*0 -1)/Pow((pi*pi/(a*a) + pi*pi/(b*b)),2))*Sin(pi*x/a)*Sin(pi*y/b)*Pow(e, -(pi*pi/(a*a) + pi*pi/(b*b))*t)";
            h1Text.Text = "0,1";
            h2Text.Text = "0,1";
            tauText.Text = "0,1";
        }

        

        private void MenuItem_Test4(object sender, RoutedEventArgs e)
        {
            fText.Text = "t*Sin(pi*x/a)*Sin(pi*y/b)";
            phiText.Text = "Sin(2*pi*x/a)*Sin(pi*y/b)";
            muText.Text = "0";
            aText.Text = "5";
            bText.Text = "5";
            solutionText.Text = "(Pow(e,(pi*pi/(a*a) + pi*pi/(b*b))*t)*((pi*pi/(a*a) + pi*pi/(b*b))*t -1)/Pow((pi*pi/(a*a) + pi*pi/(b*b)),2)-Pow(e,(pi*pi/(a*a) + pi*pi/(b*b))*0)*((pi*pi/(a*a) + pi*pi/(b*b))*0 -1)/Pow((pi*pi/(a*a) + pi*pi/(b*b)),2))*Sin(pi*x/a)*Sin(pi*y/b)*Pow(e, -(pi*pi/(a*a) + pi*pi/(b*b))*t)+Sin(2*pi*x/a)*Sin(pi*y/b)*Pow(e, -(4*pi*pi/(a*a) + pi*pi/(b*b))*t)";
            h1Text.Text = "0,1";
            h2Text.Text = "0,1";
            tauText.Text = "0,1";
        }

        private void MenuItem_Test5(object sender, RoutedEventArgs e)
        {
            fText.Text = "t*Sin(pi*x/a)*Sin(pi*y/b)+2";
            phiText.Text = "Sin(2*pi*x/a)*Sin(pi*y/b)";
            muText.Text = "2*t";
            aText.Text = "5";
            bText.Text = "5";
            solutionText.Text = "2*t + (Pow(e,(pi*pi/(a*a) + pi*pi/(b*b))*t)*((pi*pi/(a*a) + pi*pi/(b*b))*t -1)/Pow((pi*pi/(a*a) + pi*pi/(b*b)),2)-Pow(e,(pi*pi/(a*a) + pi*pi/(b*b))*0)*((pi*pi/(a*a) + pi*pi/(b*b))*0 -1)/Pow((pi*pi/(a*a) + pi*pi/(b*b)),2))*Sin(pi*x/a)*Sin(pi*y/b)*Pow(e, -(pi*pi/(a*a) + pi*pi/(b*b))*t)+Sin(2*pi*x/a)*Sin(pi*y/b)*Pow(e, -(4*pi*pi/(a*a) + pi*pi/(b*b))*t)";
            h1Text.Text = "0,1";
            h2Text.Text = "0,1";
            tauText.Text = "0,01";
        }

        private void next_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                p.nextU();
            }
            catch
            {
                MessageBox.Show("Ошибка. Проверьте начальные условия");
                return;
            }

            p2.nextU();
            //Если дано решение, рисуем
            if (solutionFlag)
            {
                showChart();
                p.calcSolution(solution);     //Считаем решение на сетке

                chart1.plotAnalytic(p);       //Рисуем
                chart1.TransformChart();

                double r = p.errorValue();    //Погрешность численного решение отностительно аналитического  
                SeriesCollection[0].Values.Add(r); //Добавляем точку на график
            }
            else
            {
                hideChart();
            }

            //Рисуем численное решение
            chart2.plotNumerical(p);
            chart2.TransformChart();

            //Считаем погрешность по рунге
            double tempR = rungeR(p, p2);
            runge = Math.Max(tempR, runge);
            SeriesCollection[1].Values.Add(tempR);

            //Обнавляем погрешность на текущий момент времени
            maxRText.Text = "R = " + p.maxR.ToString("N6");
            maxRRText.Text = "R2 = " + runge.ToString("N6");
        }

        private void Play_Click(object sender, RoutedEventArgs e)
        {
            animateFlag = !animateFlag;
            if (animateFlag)
            {
                Play.Content = "⏸";
            }
            else
            {
                Play.Content = "▶";
            }
        }

        private void MenuItem_Click(object sender, RoutedEventArgs e)
        {
            System.Windows.Application.Current.Shutdown();
        }

        private void MenuItem_Click_1(object sender, RoutedEventArgs e)
        {
            MessageBox.Show("Все функции необходимо вводить с большой буквы.\nДля возведения в степень используйте функцию Pow(a,b).\nНе забывайте ставить знак умножения.\n\nВнимание!\nПри вводе дробного числа в поле для функций испльзуйте точку.\nПри вводе дробного числа в поле для шагов сетки используйте запятую.\n\nДля разбора выражений используется NCalc. \nБолее подробно с ним можно ознакомиться в интернете.", "Правила ввода");
        }
    }
}
