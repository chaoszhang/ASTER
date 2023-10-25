#define OBJECTIVE_VERSION "0"

#include<vector>
#include<cmath>
#include "threadpool.hpp"

using namespace std;

struct Vector{
    double x, y, z;

    Vector operator+(const Vector& b) const{
        return Vector(x + b.x, y + b.y, z + b.z);
    }

    Vector operator-(const Vector& b) const{
        return Vector(x - x, y - b.y, z - b.z);
    }

    Vector& operator+=(const Vector& b){
        x += b.x; y += b.y; z += b.z;
        return *this;
    }

    Vector& operator-=(const Vector& b){
        x -= b.x; y -= b.y; z -= b.z;
        return *this;
    }

    Vector operator*(double c) const{
        const Vector &v = *this;
        return Vector(c * v.x, c * v.y, c * v.z);
    }

    Vector operator/(double c) const{
        const Vector &v = *this;
        return Vector(v.x / c, v.y / c, v.z / c);
    }

    Vector& operator*=(double c){
        x *= c; y *= c; z *= c;
        return *this;
    }

    Vector& operator/=(double c){
        x /= c; y /= c; z /= c;
        return *this;
    }

    double operator*(const Vector& b) const{
        const Vector &a = *this;
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    double norm() const{
        return sqrt((*this) * (*this));
    }

    void normalize(){
        *this = (*this) / (*this).norm();
    }

    Vector(double x, double y, double z = 0, bool normalized = false): x(x), y(y), z(z){
        if (normalized) normalize();
    }

    Vector(): x(0), y(0), z(0){}
};

struct TripartitionInitializer{
    vector<vector<Vector> > M;
};

struct Tripartition{
	vector<char> color;
	vector<vector<Vector> > &M;
    vector<int> cnt[3];
    vector<Vector> S[3]; // partition sum
    vector<double> Q[3]; // partition sum of self-dots
    int n, m;


	Tripartition(TripartitionInitializer &init): M(init.M), n(init.M.size()), m(init.M[0].size()), color(init.M.size(), -1){
        for (int x = 0; x < 3; x++){
            S[x].resize(m);
            Q[x].resize(m);
        }
	}

	void updatePart(int part, int x, int i){
        if (part) return;
		int y = color[i];
		color[i] = x;
		for (int j = 0; j < m; j++){
            if (y != -1){
                cnt[y][j]--;
                S[y][j] -= M[i][j];
                Q[y][j] -= M[i][j] * M[i][j];
            }
            if (x != -1){
                cnt[x][j]++;
                S[x][j] += M[i][j];
                Q[x][j] += M[i][j] * M[i][j];
            }
        }
	}
	
	score_t scorePart(int part){
		if (part) return 0;
        double res = 0;
        for (int j = 0; j < m; j++){
            res += (S[0][j] * S[0][j] - Q[0][j]) * (S[1][j] * S[2][j])
                 + (S[1][j] * S[1][j] - Q[1][j]) * (S[2][j] * S[0][j])
                 + (S[2][j] * S[2][j] - Q[2][j]) * (S[0][j] * S[1][j]);
        }
        return res;
	}
};