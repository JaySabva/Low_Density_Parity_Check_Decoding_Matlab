#include<bits/stdc++.h>
using namespace std;
int main()
{
    int n = 12, k = 3, u = 9;
    int H[u][n] =
    {  
    {1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0},
    {1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0},
    {0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0},
    {0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1},
    {0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1},
    {0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0},
    {1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0},
    {0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0},
    {0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1},
    };
    int DC = 0, DV = 0;
    for(int i = 0; i < u; i++)
    {
        if(H[i][0] == 1) DV++;
    }
    for(int i = 0; i < n; i++)
    {
        if(H[0][i] == 1) DC++;
    }
    int codeword[12] = {0};
    int revcode[12] = {9, 9, 9, 0, 0, 9, 0, 0, 0, 0, 0, 0};

    int CN[u][2 * DC] = {0};
    int VN[n][DV + 1] = {0};

    for (int vn = 0; vn < n; vn++)
    {
        int idx = 0;
        for(int cn = 0; cn < u; cn++)
        {
            if(H[cn][vn] == 1)
            {
                VN[vn][idx] = cn;
                idx ++;
            }
        }
    }
    for(int cn = 0; cn < u; cn++)
    {
        int idx = 0;
        for(int vn = 0; vn < n; vn++)
        {
            if(H[cn][vn] == 1)
            {
                CN[cn][idx] = vn;
                idx ++;
            }
        }
    }
    for(int vn = 0; vn < n; vn++)
    {
        VN[vn][n - 1] = revcode[vn];
    }
    for (int t = 0; t <50;t++)
    {
        for(int cn = 0; cn < u; cn++)
        {
            for(int vn = DC; vn < 2*DC; vn++)
            {
                CN[cn][vn] = VN[CN[cn][vn - DC]][n - 1];
            }
        }
    }
    int noterror = 0;
    for(int cn = 0; cn < u;cn++)
    {
        vector<int> error(CN[cn] + DC, CN[cn] + 2 * DC);
        vector<int> ers;
        for(int i = 0; i < error.size(); i++)
        {
            if(error[i] == 9)
            {
                ers.push_back(i);
            }
        }
        if(ers.size() == 0)
        {
            noterror++;
        }
        else if(ers.size() == 1)
        {
            int vnidx = CN[cn][ers[0]] - 1;
            int sum = 0;
            int era = 0;
            for (int i = DC; i < 2*DC;cn++)
            {
                if(CN[cn][i] == 1 || CN[cn][i] == 0)
                {
                    sum = sum + CN[cn][i];
                    sum = sum % 2;
                }
                else
                {
                    era = era + i;
                }
            }
            sum = sum + era;
            if(sum == 9)
            {
                VN[vnidx][n - 1] = 0;
            }
            else if(sum == 10)
            {
                VN[vnidx][n - 1] = 1;
            }
        }
        if(noterror == u)
        {
            break;
        }
    }
    int decoded[12] = {0};
    for(int i = 0; i < n;i++)
    {
        decoded[i] = VN[i][n - 1];
    }
    for(int i = 0; i < n;i++)
    {
        cout << decoded[i] << " ";
    }
}
