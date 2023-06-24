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
    vector<vector<int>> Node(u, vector<int>(n, -1));
    vector<int> decoded(n, 9);
    for (int t = 0; t < 100; t++)
    {
        if(t == 0)
        {
            for(int c = 0; c < u; c++)
            {
                for(int v = 0; v < n; v++)
                {
                    if(H[c][v] == 1)
                    {
                        Node[c][v] = revcode[v];
                    }
                }
            }
            for(int i = 0; i < n; i++)
            {
                decoded[i] = revcode[i];
            }
        }
        else
        {
            for (int c = 0; c < u; c++)
            {
                for(int v = 0; v < n; v++)
                {
                    if(Node[c][v] == 9)
                    {
                        int sum = 0;
                        for(int i = 0; i < n; i++)
                        {
                            if(i != v && Node[c][i] == 9)
                            {
                                Node[c][v] = 9;
                                goto jump;
                            }
                            else if(i != v && Node[c][i] != -1)
                            {
                                sum = sum ^ Node[c][i];
                                Node[c][v] = sum;
                            }
                        }
                    }
                }
                jump:;
            }
            for (int v = 0; v < n; v++)
            {
                for(int c = 0; c < u; c++)
                {
                    if(Node[c][v] != 9 && Node[c][v] != -1)
                    {
                        for(int i = 0; i < u; i++)
                        {
                            if(Node[i][v] == 9 && i != c)
                            {
                                Node[i][v] = Node[c][v];
                            }
                        }
                        decoded[v] = Node[c][v];
                        break;
                    }
                }
            }
        }
    }
    for(int i = 0; i < n; i++)
    {
        cout << decoded[i] << " " ;
    }
}