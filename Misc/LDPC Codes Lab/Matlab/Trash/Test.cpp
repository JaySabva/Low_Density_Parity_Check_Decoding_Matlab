#include <bits/stdc++.h>
using namespace std;
int main()
{
    int n = 9, k = 4, u = 5; // this can be derived from the H matrix
    int H[u][n] =
        {
            {1, 1, 1, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 1, 1, 1, 0, 0, 0},
            {1, 0, 0, 1, 0, 0, 1, 0, 0},
            {0, 1, 0, 1, 0, 0, 0, 1, 0},
            {0, 0, 1, 0, 0, 1, 0, 0, 1},
        };
    int dc = 0, DV = 0; // sum of row and column
    for (int i = 0; i < u; i++)
    {
        if (H[i][0] == 1)
            DV++;
    }
    for (int i = 0; i < n; i++)
    {
        if (H[0][i] == 1)
            dc++;
    }
    int revcode[9] = {9, 9, 1, 1, 9, 0, 0, 9, 9};
    vector<vector<int>> CN(u, vector<int>(2 * dc, 0));
    vector<vector<int>> VN(n, vector<int>(DV + 1, 0));

    for (int i = 0; i < n; i++)
    {
        int idx = 0;
        for (int j = 0; j < u; j++)
        {
            if (H[j][i] == 1)
            {
                VN[i][idx] = j;
                idx++;
            }
        }
    }
    for (int i = 0; i < u; i++)
    {
        int idx = 0;
        for (int j = 0; j < n; j++) 
        {
            if (H[i][j] == 1)
            {
                CN[i][idx] = j;
                idx++;
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        VN[i][n - 1] = revcode[i];
    }

    for (int t = 0; t < 50; t++)
    {
        for (int i = 0; i < n; i++)
        {
            if (VN[i][n - 1] == 9)
                ;
        }
        for (int i = 0; i < u; i++)
        {
            for (int j = dc; j < 2 * dc; j++)
            {
                CN[i][j] = VN[CN[i][j - dc]][n - 1];
            }
        }
        int error = 0;
        for (int cn = 0; cn < u; cn++)
        {
            int sumcn = accumulate(CN[cn].begin() + dc, CN[cn].begin() + 2 * dc, 0);
            vector<int> list(CN[cn].begin() + dc, CN[cn].begin() + 2 * dc);
            vector<int> erasures;
            for (int i = 0; i < list.size(); i++)
            {
                if (list[i] == 9)
                {
                    erasures.push_back(i);
                }
            }
            if (erasures.size() == 0)
            {
                error++;
            }
            else if (erasures.size() == 1)
            {
                int vnidx = CN[cn][erasures[0]];
                if((sumcn - 1) % 2 == 0)
                    VN[vnidx][n - 1] = 0;
                else
                    VN[vnidx][n - 1] = 9;
            }
        }
        if(error == u)
        {
            break;
        }
    }
            vector<int> result;
        for(int i =0 ; i < VN.size(); i++)
        {
            result.push_back(VN[i][n - 1]);
        }

    //print result
    for(int i = 0; i < n; i++)
    {
        cout << result[i] << " " ;
    }
}