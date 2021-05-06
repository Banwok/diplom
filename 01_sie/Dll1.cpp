#define BOOST_PYTHON_STATIC_LIB
#include <boost/python.hpp>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iostream>

namespace py = boost::python;

double a_p = 0.001, a_z = 0.0001, alfa_0 = 0.23, b_p = 0.01, b_z = 0.01, g_p = 0.001, q_p = 6,
q_z = 5.5, Z0 = 0, Z1 = 1, P0 = 0, P1 = 1, k_z = 0.15, k_p = 0.05, Q0 = 5, Qs = 0.2, t_porog = 125000;

size_t t = 200000;
double t0 = 0, h = 0.05, x0 = 0.2, y0_ = 0.2;
double i_py = 0, i_end = 0.012;

double fx(double Z, double P, double t, double w_q) {
  return (-(a_z + g_p * P) * Z + b_z * (Z0 - (Z0 - Z1) / (1 + exp(-(Q0 + alfa_0 * Z + Qs * sin(w_q * t) - q_z) / k_z))));
}

double fy(double Z, double P, double t, double w_q) {
  return (-a_p * P + b_p * (P0 - (P0 - P1) / (1 + exp(-(Q0 + alfa_0 * Z + Qs * sin(w_q * t) - q_p) / k_p))));
}

void init(double a_p1, double a_z1, double alfa_01, double b_p1, double b_z1, double g_p1, double q_p1, double q_z1, double Q01, double Qs1) {
  a_p = a_p1; a_z = a_z1; alfa_0 = alfa_01; 
  b_p = b_p1; b_z = b_z1; 
  g_p = g_p1; 
  q_p = q_p1;q_z = q_z1; 
  Q0 = Q01; Qs = Qs1;
}

void init_rarely(double Z01, double Z11, double P01, double P11, double k_z1, double k_p1, double t_porog1) {
  Z0 = Z01; Z1 = Z11;
  P0 = P01; P1 = P11;
  k_z = k_z1; k_p = k_p1;
  t_porog = t_porog1;
}

void init_runge(double p_x0, double p_y0, size_t p_t, double p_h, double p_time0) {
  t = p_t;
  t0 = p_time0; h = p_h;
  x0 = p_x0; y0_ = p_y0;
}

void init_wq(double i_per, double end_per) {
  i_py = i_per;
  i_end = end_per;
}

/*___________________________________________________________________________________________________________________________*/
double dot(const std::vector<double>& vec1, const std::vector<double>& vec2) {
  int size = vec1.size();
  if (size != vec2.size()) {
    std::cout << "V1 = " << vec1.size();
    std::cout << "\nV2 = " << vec2.size();
    throw 1;
  }
  double sum = 0;

  for (int i = 0; i < size; ++i)
    sum += vec1[i] * vec2[i];

  return sum;
}

std::vector<double> addit_vec(std::vector<double>& vec1, std::vector<double>& vec2, int fl) {
  int size = vec1.size();
  if (size != vec2.size()) {
    std::cout << "V1 = " << vec1.size();
    std::cout << "\nV2 = " << vec2.size();
    throw 1;
  }

  std::vector<double> tmp = vec2;
  if (fl == 1) {
    for (size_t i = 0; i < size; ++i)
      tmp[i] *= 0.5;
  }

  auto res = std::vector<double>(size);
  std::transform(vec1.begin(), vec1.end(), tmp.begin(), res.begin(), std::plus<double>());
  return res;
}

std::vector<double> minus_vec(std::vector<double>& vec1, std::vector<double>& vec2, double fl) {
  int size = vec1.size();
  if (size != vec2.size())
    throw 1;

  std::vector<double> tmp = vec2;
  for (size_t i = 0; i < size; ++i)
    tmp[i] *= fl;

  auto res = std::vector<double>(size);
  std::transform(vec1.begin(), vec1.end(), tmp.begin(), res.begin(), std::minus<double>());
  return res;
}

double norm(std::vector<double>& vec) {
  double a = vec[0];
  double b = vec[1];
  double res = sqrt(pow(a, 2) + pow(b, 2));
  return res;
}

std::vector<double> vec_Jac_X(const std::vector<double>& dfx, double Z, double P, double t, double w_q) {
  double fx_z_coff = (exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_z) / k_z) + 1);
  double fy_z_coff = (exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_p) / k_p) + 1);

  double fx_Z = -P * g_p - a_z - alfa_0 * b_z * (Z0 - Z1) * exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_z) / k_z) / (k_z * fx_z_coff * fx_z_coff);
  double fx_P = -Z * g_p;
  double fy_Z = -alfa_0 * b_p * (P0 - P1) * exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_p) / k_p) / (k_p * fy_z_coff * fy_z_coff);
  double fy_P = -a_p;

  std::vector<double> d_fx{ fx_Z , fx_P };
  std::vector<double> d_fy{ fy_Z, fy_P };

  double d1_X = dot(d_fx, dfx) * h;
  double d2_X = dot(d_fy, dfx) * h;
  std::vector<double> res = { d1_X, d2_X };
  return res;
}

std::vector<double> vec_Jac_Y(std::vector<double>& dfy, double Z, double P, double t, double w_q) {
  double fx_z_coff = (exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_z) / k_z) + 1);
  double fy_z_coff = (exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_p) / k_p) + 1);

  double fx_Z = -P * g_p - a_z - alfa_0 * b_z * (Z0 - Z1) * exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_z) / k_z) / (k_z * fx_z_coff * fx_z_coff);
  double fx_P = -Z * g_p;
  double fy_Z = -alfa_0 * b_p * (P0 - P1) * exp((-Q0 - Qs * sin(t * w_q) - Z * alfa_0 + q_p) / k_p) / (k_p * fy_z_coff * fy_z_coff);
  double fy_P = -a_p;

  std::vector<double> d_fx{ fx_Z , fx_P };
  std::vector<double> d_fy{ fy_Z, fy_P };

  double d1_Y = dot(d_fx, dfy) * h;
  double d2_Y = dot(d_fy, dfy) * h;
  std::vector<double> res = { d1_Y, d2_Y };
  return res;
}

void RK_lyapunov(std::vector<double>& a, std::vector<double>& b, std::vector<double>& lya2, std::vector<double>& lya1,  double x0, double y0,
  size_t t, double h, double time0, double w_q, size_t& iter, size_t &j) {
  size_t size = static_cast<int>(t / h + 1);
  std::vector<double> x(size);
  std::vector<double> y(size);

  double yt = 0;
  double xt = 0;
  double l1 = 0;
  double l2 = 0;

  double t0 = time0;
  size_t i = 0;
  x[i] = x0;
  y[i] = y0;

  std::vector<double> dfx{ 1, 0 };
  std::vector<double> dfy{ 0, 1 };
  std::vector<double> ort1(2);
  std::vector<double> ort2(2);

  for (t0; t0 < t; t0 += h) {
    double k1_y =  fy(x[i], y[i], t0, w_q);
    double k2_y = fy(x[i] + k1_y / 2.0, y[i] + k1_y / 2.0, t0, w_q);
    double k3_y = fy(x[i] + k2_y / 2.0, y[i] + k2_y / 2.0, t0, w_q);
    double k4_y = fy(x[i] + k3_y, y[i] + k3_y, t0, w_q);
    yt = y[i] + h *(k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
    y[i + 1] = yt;

    double k1_x = fx(x[i], y[i], t0, w_q);
    double k2_x = fx(x[i] + k1_x / 2.0, y[i] + k1_x / 2.0, t0, w_q);
    double k3_x = fx(x[i] + k2_x / 2.0, y[i] + k2_x / 2.0, t0, w_q);
    double k4_x =  fx(x[i] + k3_x, y[i] + k3_x, t0, w_q);
    xt = x[i] + h * (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
    x[i + 1] = xt;

    std::vector<double> k11_x = vec_Jac_X(dfx, x[i], y[i], t0, w_q);
    std::vector<double> tmp = addit_vec(dfx, k11_x, 1);
    std::vector<double> k22_x = vec_Jac_X(tmp, x[i] + k1_x / 2.0 * h, y[i] + k1_x / 2.0 * h, t0, w_q);
    tmp = addit_vec(dfx, k22_x, 1);
    std::vector<double> k33_x = vec_Jac_X(tmp, x[i] + k2_x / 2.0 * h, y[i] + k2_x / 2.0 * h, t0, w_q);
    tmp = addit_vec(dfx, k33_x, 0);
    std::vector<double> k44_x = vec_Jac_X(tmp, x[i] + k3_x * h, y[i] + k3_x * h, t0, w_q);

    std::vector<double> k11_y = vec_Jac_Y(dfy, x[i], y[i], t0, w_q);
    std::vector<double> tmpy = addit_vec(dfy, k11_y, 1);
    std::vector<double> k22_y = vec_Jac_Y(tmpy, x[i] + k1_y / 2.0 * h, y[i] + k1_y / 2.0 * h, t0, w_q);
    tmpy = addit_vec(dfy, k22_y, 1);
    std::vector<double> k33_y = vec_Jac_Y(tmpy, x[i] + k2_y / 2.0 * h, y[i] + k2_y / 2.0 * h, t0, w_q);
    tmpy = addit_vec(dfy, k33_y, 0);
    std::vector<double> k44_y = vec_Jac_Y(tmpy, x[i] + k3_y * h, y[i] + k3_y * h, t0, w_q);

    dfx = { dfx[0] + (k11_x[0] + 2.0 * k22_x[0] + 2.0 * k33_x[0] + k44_x[0]) / 6.0, dfx[1] + (k11_x[1] + 2.0 * k22_x[1] + 2.0 * k33_x[1] + k44_x[1]) / 6.0 };
    dfy = { dfy[0] + (k11_y[0] + 2.0 * k22_y[0] + 2.0 * k33_y[0] + k44_y[0]) / 6.0, dfy[1] + (k11_y[1] + 2.0 * k22_y[1] + 2.0 * k33_y[1] + k44_y[1]) / 6.0};

    ort1 = dfx;
    double n1 = norm(ort1);
    l1 += log(n1);
    dfx = { ort1[0] / n1, ort1[1] / n1 };
    
    double dot_ = dot(dfy, dfx);
    ort2 = minus_vec(dfy, dfx, dot_);
    double n2 = norm(ort2);
    l2 += log(n2);
    dfy = { ort2[0] / n2, ort2[1] / n2 };

    if ((t0 > t_porog) && ((x[i - 2] - x[i - 1]) * (x[i - 1] - x[i]) < 0) && ((x[i - 2] - x[i - 1]) < 0)) {
      a[iter] = x[i - 1];
      b[iter] = w_q;
      ++iter;
    }
    ++i;
  }
  double lk1, lk2;
  double TT = static_cast<double>(t);
  lk1 = (l1 / TT) / log(2);
  lk2 = (l2 / TT) / log(2);
  lya1[j] = lk1;
  lya2[j] = lk2;
}

void lyapunov_solution(py::list localM, py::list w_s, py::list ly1, py::list ly2, py::list w) {
  size_t iter = 0;
  size_t j = 0;
  std::vector<double> lmax(9000000);
  std::vector<double> w_SS(9000000);
  int ili = static_cast<int>(i_py * 100000);
  int Iend = static_cast<int>(i_end * 100000);
  std::vector<double> lyapun1(Iend - ili);
  std::vector<double> lyapun2(Iend - ili);
  std::vector<double> wwwq(Iend - ili);

  double start = omp_get_wtime();
#pragma omp parallel for shared(iter, j)
  for (int i = ili; i < Iend; i += 1) {
    double w_q = static_cast<double>(i) / 100000;
    RK_lyapunov(lmax, w_SS, lyapun1, lyapun2, x0, y0_, t, h, t0, w_q, iter, j);
    wwwq[j] = w_q;
    ++j;
  }
  double end = omp_get_wtime();
  printf("Parallel lyapunov work time %f seconds\n", end - start);

  lmax.erase(lmax.begin() + iter, lmax.end());
  w_SS.erase(w_SS.begin() + iter, w_SS.end());

  std::vector<double>(lmax).swap(lmax);
  std::vector<double>(w_SS).swap(w_SS);

  lmax.erase(std::remove(lmax.begin(), lmax.end(), 0.0), lmax.end());
  w_SS.erase(std::remove(w_SS.begin(), w_SS.end(), 0.0), w_SS.end());

  for (size_t k = 0; k < lmax.size(); ++k)
    localM.append(lmax[k]);

  for (size_t j = 0; j < w_SS.size(); ++j)
    w_s.append(w_SS[j]);

  for (size_t l_i = 0; l_i < lyapun1.size(); ++l_i)
    ly1.append(lyapun1[l_i]);

  for (size_t l_j = 0; l_j < lyapun2.size(); ++l_j)
    ly2.append(lyapun2[l_j]);

  for (size_t j = 0; j < wwwq.size(); ++j)
    w.append(wwwq[j]);
}

/*___________________________________________________________________________________________________________________________*/

void RKutt(std::vector<double> &a, std::vector<double>&b,double x0, double y0, size_t t, double h, double time0, double w_q, size_t &iter) {
  size_t size = static_cast<int>(t / h + 1);
  std::vector<double> x(size);
  std::vector<double> y(size);
  double yt = 0;
  double xt = 0;
  double t0 = time0;
  size_t i = 0;
  x[i] = x0;
  y[i] = y0;
  
  for (t0; t0 < t; t0 += h) {
    double k1_y = fy(x[i], y[i], t0, w_q);
    double k2_y = fy(x[i] + k1_y / 2.0, y[i] + k1_y / 2.0, t0, w_q);
    double k3_y = fy(x[i] + k2_y / 2.0, y[i] + k2_y / 2.0, t0, w_q);
    double k4_y = fy(x[i] + k3_y, y[i] + k3_y, t0, w_q);
    yt = y[i] + h * (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
    y[i + 1] = yt;

    double k1_x = fx(x[i], y[i], t0, w_q);
    double k2_x = fx(x[i] + k1_x / 2.0, y[i] + k1_x / 2.0, t0, w_q);
    double k3_x = fx(x[i] + k2_x / 2.0, y[i] + k2_x / 2.0, t0, w_q);
    double k4_x = fx(x[i] + k3_x, y[i] + k3_x, t0, w_q);
    xt = x[i] + h * (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
    x[i + 1] = xt;

    if ((t0 > t_porog) && ((x[i - 2] - x[i - 1]) * (x[i - 1] - x[i]) < 0) && ((x[i - 2] - x[i - 1]) < 0)) {
      a[iter] = x[i - 1];
      b[iter] = w_q;
      ++iter;
    }
    ++i;
  }
}

void bMotherFurkation(py::list localM, py::list w_s) {
  size_t iter = 0;
  std::vector<double> lmax(9000000);
  std::vector<double> w_SS(9000000);
  int ili = static_cast<int>(i_py * 100000);
  
  int Iend = static_cast<int>(i_end * 100000);
  double start = omp_get_wtime();
#pragma omp parallel for shared(iter)
  for (int i = ili; i < Iend; i += 1) {
    double w_q = static_cast<double>(i) / 100000;
    RKutt(lmax, w_SS, x0, y0_, t, h, t0, w_q, iter);
  }
  double end = omp_get_wtime();
  printf("Parallel work time %f seconds\n", end - start);

  lmax.erase(lmax.begin() + iter, lmax.end());
  w_SS.erase(w_SS.begin() + iter, w_SS.end());
  
  std::vector<double>(lmax).swap(lmax);
  std::vector<double>(w_SS).swap(w_SS);

  lmax.erase(std::remove(lmax.begin(), lmax.end(), 0.0), lmax.end());
  w_SS.erase(std::remove(w_SS.begin(), w_SS.end(), 0.0), w_SS.end());

  for (size_t k = 0; k < lmax.size(); ++k)
    localM.append(lmax[k]);

  for (size_t j = 0; j < w_SS.size(); ++j)
    w_s.append(w_SS[j]);
}

/* ____________________________________________________________________________________________*/

void RK_wiki(std::vector<double>& x, std::vector<double>& y, std::vector<double> &time,double x0, double y0, 
                size_t t, double h, double time0, double w_q) {
  double yt = 0;
  double xt = 0;
  double t0 = time0;
  size_t i = 0;
  time[i] = t0;
  x[i] = x0;
  y[i] = y0;

  for (t0; t0 < t; t0 += h) {
    double k1_y = fy(x[i], y[i], t0, w_q);
    double k2_y = fy(x[i] + k1_y / 2.0, y[i] + k1_y / 2.0, t0, w_q);
    double k3_y = fy(x[i] + k2_y / 2.0, y[i] + k2_y / 2.0, t0, w_q);
    double k4_y = fy(x[i] + k3_y, y[i] + k3_y, t0, w_q);
    yt = y[i] + h * (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
    y[i + 1] = yt;

    double k1_x = fx(x[i], y[i], t0, w_q);
    double k2_x = fx(x[i] + k1_x / 2.0, y[i] + k1_x / 2.0, t0, w_q);
    double k3_x = fx(x[i] + k2_x / 2.0, y[i] + k2_x / 2.0, t0, w_q);
    double k4_x = fx(x[i] + k3_x, y[i] + k3_x, t0, w_q);
    xt = x[i] + h * (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
    x[i + 1] = xt;
    time[i + 1] = t0;
    ++i;
  }
}

void retuRelease(py::list funX, py::list funY, py::list funT, double w_q) {
  size_t size = static_cast<int>(t / h + 1);
  std::vector<double> x(size);
  std::vector<double> y(size);
  std::vector<double> time(size);

  RK_wiki(x, y, time, x0, y0_, t, h, t0, w_q);

  for (size_t i = 0; i < x.size(); ++i)
    funX.append(x[i]);

  for (size_t i = 0; i < y.size(); ++i)
    funY.append(y[i]);

  for (size_t i = 0; i < time.size(); ++i)
    funT.append(time[i]);
}

/*_____________________________________________________________________________________________________*/

void RKutt_seq(std::vector<double>& a, std::vector<double>& b, double x0, double y0,
                size_t t, double h, double time0, double w_q, size_t& iter, double& retX, double& retY) {
  size_t size = static_cast<int>(t / h + 1);
  std::vector<double> x(size);
  std::vector<double> y(size);
  std::vector<double> locakMaxim;
  double yt = 0;
  double xt = 0;
  double t0 = time0;
  size_t i = 0;
  x[i] = x0;
  y[i] = y0;

  for (t0; t0 < t; t0 += h) {
    double k1_y = fy(x[i], y[i], t0, w_q);
    double k2_y = fy(x[i] + k1_y / 2.0, y[i] + k1_y / 2.0, t0, w_q);
    double k3_y = fy(x[i] + k2_y / 2.0, y[i] + k2_y / 2.0, t0, w_q);
    double k4_y = fy(x[i] + k3_y, y[i] + k3_y, t0, w_q);
    yt = y[i] + h * (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;
    y[i + 1] = yt;

    double k1_x = fx(x[i], y[i], t0, w_q);
    double k2_x = fx(x[i] + k1_x / 2.0, y[i] + k1_x / 2.0, t0, w_q);
    double k3_x = fx(x[i] + k2_x / 2.0, y[i] + k2_x / 2.0, t0, w_q);
    double k4_x = fx(x[i] + k3_x, y[i] + k3_x, t0, w_q);
    xt = x[i] + h * (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
    x[i + 1] = xt;
    ++i;
  }
  retX = xt;
  retY = yt;
}

void seq_Furkation(py::list localM, py::list w_s) {
  size_t iter = 0;
  std::vector<double> lmax(9000000);
  std::vector<double> w_SS(9000000);
  int ili = static_cast<int>(i_py * 100000);
  int Iend = static_cast<int>(i_end * 100000);
  double rX = x0;
  double rY = y0_;

  double start = omp_get_wtime();
  for (int i = ili; i < Iend; i += 1) {
    double w_q = static_cast<double>(i) / 100000;
    RKutt_seq(lmax, w_SS, x0, y0_, t, h, t0, w_q, iter, rX, rY);
  }
  double end = omp_get_wtime();
  printf("Sequential work time %f seconds\n", end - start);

  lmax.erase(lmax.begin() + iter, lmax.end());
  w_SS.erase(w_SS.begin() + iter, w_SS.end());

  std::vector<double>(lmax).swap(lmax);
  std::vector<double>(w_SS).swap(w_SS);

  lmax.erase(std::remove(lmax.begin(), lmax.end(), 0.0), lmax.end());
  w_SS.erase(std::remove(w_SS.begin(), w_SS.end(), 0.0), w_SS.end());

  for (size_t k = 0; k < lmax.size(); ++k)
    localM.append(lmax[k]);

  for (size_t j = 0; j < w_SS.size(); ++j)
    w_s.append(w_SS[j]);
}


BOOST_PYTHON_MODULE(Dll1)
{
  py::def("init", init);
  py::def("init_rarely", init_rarely);
  py::def("init_runge", init_runge);
  py::def("init_wq", init_wq);

  py::def("retuRelease", retuRelease);

  py::def("seq_Furkation", seq_Furkation);
  py::def("bMotherFurkation", bMotherFurkation);
  py::def("lyapunov_solution", lyapunov_solution);

}
