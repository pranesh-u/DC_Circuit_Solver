#include<bits/stdc++.h>
#include"eigen-3.4.0/Eigen/Dense"
#include<fstream>
using namespace std;
using namespace Eigen;

typedef struct ComponentParameters{
    string name;        // Name of Component
    int node1,node2;    // Starting and ending nodes
    double mag;         // Value of component
    double Y;           // Admittance of component
    double voltage;     // Value of Voltage Source
    double current;     // Value of Current Source
    int branch_no;      // Branch Number

}Component;

ifstream input;
ofstream output;
int VS_num = 0,CS_num = 0,R_num = 0;
vector<Component> components;

//Reads the given netlist and adds it to the components array;
void component_read()
{
    string compo_type;
    input.open("netlist.txt");
    if(input.is_open())
    {
        while(!input.eof())
        {
            input >> compo_type;
            //If the component is a voltage source
            if(compo_type[0] == 'V')
            {
                string n1, n2, voltage; 
                input >> n1 >> n2 >> voltage;
                Component voltage_source;
                voltage_source.name = compo_type;
                voltage_source.node1 = stoi(n1);
                voltage_source.node2 = stoi(n2);
                voltage_source.voltage = stod(voltage);
                voltage_source.mag = 1;
                // For Voltage Source, if Y = infinity, then the answer is NaN
                voltage_source.Y = 1;   
                components.push_back(voltage_source);
                VS_num++;
            }
            //If the component is a current source
            else if(compo_type[0] == 'I')
            {
                string n1, n2, current; 
                input >> n1 >> n2 >> current;
                Component current_source;
                current_source.name = compo_type;
                current_source.node1 = stoi(n1);
                current_source.node2 = stoi(n2);
                current_source.current = stod(current);
                current_source.mag = numeric_limits<double>::infinity();
                current_source.Y = pow(current_source.mag,-1);
                components.push_back(current_source);
                CS_num++;
            }
            //If the component is a resistor
            else if(compo_type[0] == 'R')
            {
                string n1, n2, mag; 
                input >> n1 >> n2 >> mag;
                Component resistance;
                resistance.name = compo_type;
                resistance.node1 = stoi(n1);
                resistance.node2 = stoi(n2);
                resistance.mag = stod(mag);
                resistance.Y = pow(resistance.mag,-1);
                components.push_back(resistance);
                R_num++;
            }
        }
    }
    else
    {
        cout << "File doesn't exist" << endl;
        return;
    }
}

int no_of_nodes(vector<Component> &components)
{
    int nodes;
    for(int i = 0;i < components.size();i++)
    {
        if(components[i].node1 > nodes)
            nodes = components[i].node1;
        if(components[i].node2 > nodes)
            nodes = components[i].node2;
    }
    return nodes;
}

int main()
{
    component_read();
    vector<Component> components_sorted;
    int nodes = no_of_nodes(components);
    int BranchesCount = components.size();
    for(int i = 0;i < components.size();i++)
    {
        components[i].branch_no = i+1;
    }
    // Declaration and initalisation of required matrices
    Eigen::MatrixXd Im(nodes,BranchesCount);                 // Reduced Incidence Matrix
    Eigen::MatrixXd Ys(BranchesCount,BranchesCount);         // Admittance Matrix
    Eigen::MatrixXd Vs(BranchesCount,1);                     // Voltage Source Matrix
    Eigen::MatrixXd Is(BranchesCount,1);                     // Current Source Matrix
    Eigen::MatrixXd Vn(nodes,1);                             // Node Voltages Matrix
    for(int i = 0;i < nodes;i++)
    {
        for(int j = 0;j < BranchesCount;j++)
        {
            Im(i,j) = 0;
        }
    }
    for(int i = 0;i < BranchesCount;i++)
    {
        for(int j = 0;j < BranchesCount;j++)
        {
            Ys(i,j) = 0;
        }
        Vs(i,0) = 0;
        Is(i,0) = 0;
    }

    //Determining Voltage Source Matrix and Current Source Matrix
    for(int i = 0;i < components.size();i++)
    {
        if(components[i].name[0] == 'V')
        {
            Vs(components[i].branch_no-1,0) = -components[i].voltage;
        }
        if(components[i].name[0] == 'I')
        {
            Is(components[i].branch_no-1,0) = -components[i].current;
        }
    }
    
    // Determining Admittance Matrix
    for(int i = 0;i < components.size();i++)
    {
        Ys(components[i].branch_no-1,components[i].branch_no-1) = components[i].Y;
    }

    //Determining Incidence Matrix
    for(int i = 0;i < components.size();i++)
    {
            if(components[i].node1 == 0)
            {
                Im(components[i].node2-1,components[i].branch_no-1) = -1;
            }
            else if(components[i].node2 == 0)
            {
                Im(components[i].node1-1,components[i].branch_no-1) = 1;
            }
            else
            {
                Im(components[i].node1-1,components[i].branch_no-1) = 1;
                Im(components[i].node2-1,components[i].branch_no-1) = -1;
            }
    }

    cout << "Voltage Source Matrix : " << endl;
    cout << Vs << endl;
    cout << "Current Source Matrix : " << endl;
    cout << Is << endl;
    cout << "Admittance Matrix : " << endl;
    cout << Ys << endl;
    cout << "Incidence Matrix : " << endl;
    cout << Im << endl;

    //Master Equation is A*Y*At*Vn = A(Y*V - I)
    Eigen::MatrixXd LHS(nodes,nodes);
    LHS = Im*Ys*(Im.transpose());
    Eigen::MatrixXd RHS(nodes,1);
    RHS = Im*(Ys*Vs - Is);

    Vn = LHS.colPivHouseholderQr().solve(RHS);
    cout << "Node Voltages are : " << endl;
    cout << Vn << endl;
    return 0;
}