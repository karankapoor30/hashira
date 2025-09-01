#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iomanip>


long long baseToInt(const std::string& s, int base) {
    long long res = 0;
    for (char c : s) {
        int d;
        if ('0' <= c && c <= '9') d = c - '0';
        else if ('a' <= c && c <= 'z') d = c - 'a' + 10;
        else if ('A' <= c && c <= 'Z') d = c - 'A' + 10;
        else d = 0;
        res = res * base + d;
    }
    return res;
}


std::vector<double> solvePoly(const std::vector<double>& roots, int degree) {
    int n = degree + 1;
    std::vector<std::vector<double>> A(n, std::vector<double>(n + 1, 0.0));
    for (int i = 0; i < n; ++i) {
        double x = roots[i];
        for (int j = 0; j < n; ++j) {
            A[i][j] = pow(x, degree - j);
        }

        A[i][n] = 0.0;
    }

    A[0][0] = 1.0;
    A[0][n] = -pow(roots[0], degree);

    for (int i = 0; i < n; ++i) {

        int pivot = i;
        for (int j = i + 1; j < n; ++j)
            if (fabs(A[j][i]) > fabs(A[pivot][i]))
                pivot = j;
        if (pivot != i)
            std::swap(A[i], A[pivot]);

        double fac = A[i][i];
        if (fabs(fac) < 1e-12) continue;
        for (int j = i; j <= n; ++j)
            A[i][j] /= fac;

        for (int j = i + 1; j < n; ++j) {
            double mul = A[j][i];
            for (int k = i; k <= n; ++k)
                A[j][k] -= mul * A[i][k];
        }
    }

    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = A[i][n];
        for (int j = i + 1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
    }
    return x;
}


std::vector<long long> parseRoots(const std::string& json, int& k) {
    std::vector<long long> roots;
    size_t pos = 0;
    int n = 0;
    while ((pos = json.find("\"n\"", pos)) != std::string::npos) {
        size_t npos = json.find(":", pos);
        n = std::stoi(json.substr(npos + 1, json.find(",", npos) - npos));
        pos = npos;
        break;
    }
    pos = 0;
    while ((pos = json.find("\"k\"", pos)) != std::string::npos) {
        size_t kpos = json.find(":", pos);
        k = std::stoi(json.substr(kpos + 1, json.find(",", kpos) - kpos));
        pos = kpos;
        break;
    }
    for (int i = 1; i <= n; ++i) {
        std::string idx = "\"" + std::to_string(i) + "\"";
        size_t ipos = json.find(idx);
        if (ipos == std::string::npos) continue;
        size_t bpos = json.find("\"base\"", ipos);
        size_t vpos = json.find("\"value\"", ipos);
        size_t cpos = json.find(":", bpos);
        size_t cb = json.find(",", cpos);
        std::string sbase = json.substr(cpos + 1, cb - cpos - 1);
        int base = std::stoi(sbase);
        cpos = json.find(":", vpos);
        cb = json.find("\n", cpos);
        std::string sval = json.substr(cpos + 1, cb - cpos - 1);
        sval.erase(std::remove(sval.begin(), sval.end(), '\"'), sval.end());
        sval.erase(std::remove(sval.begin(), sval.end(), ','), sval.end());
        sval.erase(std::remove_if(sval.begin(), sval.end(), ::isspace), sval.end());
        roots.push_back(baseToInt(sval, base));
    }
    return roots;
}

int main() {
    std::ostringstream buf;
    buf << std::cin.rdbuf();
    std::string input = buf.str();

    int k = 0;
    std::vector<long long> roots = parseRoots(input, k);

    int degree = k - 1;
    if (roots.size() < k) {
        std::cerr << "Not enough roots provided to solve for degree " << degree << std::endl;
        return 1;
    }

    std::vector<double> rootsd;
    for (int i = 0; i < k; ++i)
        rootsd.push_back((double)roots[i]);

    std::vector<double> coeffs = solvePoly(rootsd, degree);

    
    std::cout << std::setprecision(8) << std::fixed;
    std::cout << "Polynomial coefficients: ";
    for (int i = 0; i < coeffs.size(); ++i) {
        std::cout << coeffs[i];
        if (i + 1 < coeffs.size()) std::cout << ", ";
    }
    std::cout << std::endl;
    return 0;
}
