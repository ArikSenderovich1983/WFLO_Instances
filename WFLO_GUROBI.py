import sys, os, datetime
import numpy as np
#from dwave_qbsolv import QBSolv
import neal
import itertools
import math
import pymzn
from gurobipy import *
#from ortools.sat.python import cp_model
import time
from pyqubo import Binary, Array
import matplotlib.pyplot as plt
import pandas as pd
import sys
#import psutil

# calculating absolute distance between two points
def abs_distance(x, y):
    return abs(x - y)


def euclidean_distance(x1, x2, y1, y2):
    return pow(pow(x1 - x2, 2) + pow(y1 - y2, 2), 0.5)


# Checking and setting wake influence between pairs of cells
# This function is modified by Jason
def check_wake(p_i, p_j):
    alpha = 0.09437
    rotor_radius = 27.881

    #xdistance = abs_distance(p_i[0], p_j[0])
    if p_j[1] > p_i[1]:
        ydistance = abs_distance(p_i[1], p_j[1])

        radius = rotor_radius + (alpha * ydistance)

        xmin = p_j[0] - radius
        xmax = p_j[0] + radius
    else:
        return False
        
    if (xmin < p_i[0] and xmax > p_i[0]):
        return True
    else:
        return False


# Calculating velocity factor between two cells due to wakes
def calc_velocity(loc_y_i, loc_y_j, u_inf):
    alpha = 0.09437
    a = 0.326795
    rotor_radius = 27.881

    # todo: check proper calculation of wind
    ydistance = abs_distance(loc_y_i,
                             loc_y_j)  # euclidean_distance(0,0,loc_y_i, loc_y_j) #abs_distance(loc_y_i, loc_y_j)
    denom = pow((alpha * ydistance / rotor_radius) + 1, 2)

    return u_inf * (1 - (2 * a / denom))


def is_in_grid(x, y, max_x, max_y):
    if (x < 0 or x > max_x) or (y < 0 or y > max_y):
        return False
    else:
        return True



# setting cell locations:
def calc_location(axis_, step):
    location_x = []
    location_y = []
    for x in range(0, axis_):
        for y in range(0, axis_):
            location_x.append(x * step)
            location_y.append(y * step)
    return location_x, location_y

def set_wind_regime(n_wind):
    wind_speeds = []
    prob_l = []
    if n_wind == 36:
        wind_speeds = np.array([8] * n_wind)
        wind_speeds = np.append(wind_speeds, [12] * n_wind)
        wind_speeds = np.append(wind_speeds, [17] * n_wind)

        prob_l = np.array([float(0.0047)] * n_wind)
        prob_l = np.append(prob_l, [0.008] * (n_wind - 9))  # 9 special regimes
        prob_l = np.append(prob_l, [0.01, 0.012, 0.0145, 0.014, 0.019, 0.014, 0.0145, 0.012, 0.01])
        prob_l = np.append(prob_l, [0.011] * (n_wind - 9))
        prob_l = np.append(prob_l, [0.013, 0.017, 0.0185, 0.03, 0.035, 0.03, 0.0185, 0.017, 0.013])
    elif n_wind ==1:


        wind_speeds = np.array([12] * n_wind)

        prob_l = np.array([1])

    return wind_speeds, prob_l

def calc_coefficients(location_x, location_y, wind_speeds, prob_l, n):
    # Epsilon contains pairs of cells where Turbines cannot be placed due to proximity.
    # Smaller than dist x step (meters) is not allowed. Absolute distance is used.
    # Larger Epsilon is allowed to constrain the problem more -- Jason
    Epsilon = []
    #dist = 1  # 100 meters between turbines (in some papers it is 200 meters so dist should be 2)
    for i in range(n):
        for j in range(n):
            if ((abs_distance(location_x[i], location_x[j]) < 200) and (
                    abs_distance(location_y[i], location_y[j]) < 200)) and (i != j) and (j, i) not in Epsilon:
                Epsilon.append((i, j))
    len(Epsilon)

    max_x = max(location_x)
    max_y = max(location_y)
    velocity = np.zeros((n, n, len(wind_speeds)))
    energy_coef = np.zeros((n, n, len(wind_speeds)))
    ss_coef = np.zeros((n, n, len(wind_speeds)))
    U = {}
    # per each cell, we want to consider the locations that if turbine is placed there it will result
    # in a wake
    for l in range(0, len(wind_speeds)):
        for i in range(n):
            U[(i, l)] = []
            for j in range(n):
                # print(str((i,j)))
                # default influence is zero (no wake)
                energy_coef[i, j, l] = 0
                if i == j:
                    continue
                # assuming wind is going south to north, wakes will happen only on the y axis bottom up
                # if location_y[i]>=location_y[j]:
                # if true, then i is northern to j so it is influenced by it but does not influence it (j is upstream to i)
                #    continue
                # else:
                # potential wake
                temp_x_i = location_x[i] * math.cos(((l % n_wind) * angle) * math.pi / 180) - location_y[i] * math.sin(
                    ((l % n_wind) * angle) * math.pi / 180)

                temp_x_j = location_x[j] * math.cos(((l % n_wind) * angle) * math.pi / 180) - location_y[j] * math.sin(
                    ((l % n_wind) * angle) * math.pi / 180)

                temp_y_i = location_x[i] * math.sin(((l % n_wind) * angle) * math.pi / 180) + location_y[i] * math.cos(
                    ((l % n_wind) * angle) * math.pi / 180)
                temp_y_j = location_x[j] * math.sin(((l % n_wind) * angle) * math.pi / 180) + location_y[j] * math.cos(
                    ((l % n_wind) * angle) * math.pi / 180)

                # print(temp_x_i, temp_x_j, temp_y_i, temp_y_j)
                # diff = euclidean_distance(location_x[i], location_x[j], location_y[i], location_y[j]) -\
                #             euclidean_distance(temp_x_i, temp_x_j, temp_y_i, temp_y_j)

                if True:  # is_in_grid(temp_x_i, temp_y_i, max_x, max_y) and is_in_grid(temp_x_j, temp_y_j, max_x, max_y):
                    # print((l,i,j))

                    prob = prob_l[l]
                    u_inf = wind_speeds[l]
                    if check_wake((temp_x_i, temp_y_i), (temp_x_j, temp_y_j)):
                        # if check_wake((0,0), (temp_x_j, temp_y_j)):
                        # if there is a wake we calculate the velocity impact and the energy loss
                        U[(i, l)].append(j)
                        velocity[i, j, l] = calc_velocity(temp_y_i, temp_y_j, u_inf)
                        energy_coef[i, j, l] = 0.33 * (pow(u_inf, 3) - pow(velocity[i, j, l], 3))
                        ss_coef[i, j, l] = pow((1 - velocity[i, j, l] / wind_speeds[l]), 2)

    aggregated_coef = np.zeros((n, n))
    # prob_l = np.array([float(1)]*n_wind)
    # prob_l = np.true_divide(prob_l,n_wind)

    for i in range(n):
        for j in range(n):
            for l in range(0, len(wind_speeds)):
                if j in U[(i, l)]:
                    aggregated_coef[i, j] += prob_l[l] * energy_coef[i, j, l]

    return aggregated_coef, ss_coef, Epsilon, U

# In[8]:

def calc_energy_full(n, wind_speeds, prob_l, X, ss_coef):
    Energy = {}
    N = range(n)
    for i in N:
        for l in range(len(wind_speeds)):
            # only wake!!!
            Energy[(i, l)] = 0.33 * prob_l[l] * X[i]
            temp = 0
            for j in N:
                if j in U[(i, l)]:
                    temp += X[j] * ss_coef[i, j, l]  # pow(1 - velocity[i, j, l] / wind_speeds[l], 2)
            Energy[(i, l)] *= pow(wind_speeds[l] * (1 - pow(temp, 0.5)), 3)
    total_energy = sum([Energy[(i, l)] for i in N for l in range(len(wind_speeds))])
    return total_energy
### Gurobi version here.
# Benchmarking against MILP (Gurobi)


# Determine what to output!
def call_back(model, where):
    '''
    if where == GRB.Callback.MIP:
        runt = int(model.cbGet(GRB.Callback.RUNTIME))
        nodecnt = model.cbGet(GRB.Callback.MIP_NODCNT)
        solcnt = model.cbGet(GRB.Callback.MIP_SOLCNT)
        objbst = model.cbGet(GRB.Callback.MIP_OBJBST)
        objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
        objgap = abs(objbst - objbnd) / (0.0000000001 + abs(objbnd))
        intermediate_.append((runt, objbst, objbnd, objgap, solcnt, nodecnt))
    '''
    if where == GRB.Callback.MIPSOL:
        try:
            x = model.cbGetSolution(model.getVars()[0:n])
            # make a list of edges selected in the solution
            nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
            obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)
            solcnt = model.cbGet(GRB.Callback.MIPSOL_SOLCNT)
            runt = int(model.cbGet(GRB.Callback.RUNTIME))
            print('Runtime is: ' + str(runt))
            objbnd = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
            objbst = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
            objgap = abs(objbst - objbnd) / (0.0000000001 + abs(objbnd))
            # I have appended nodecnt as another attribute -- Jason
            # For standard instances, no optimality gap condition of termination
            # For large instances, terminate when the optimality gap is < 0.1
            #if objgap < 0.0000000001:
                #print('Stop early - 10% gap achieved')
                #model.terminate()
            print('**** New solution at node %d, obj %g, sol %d, ' \
            '****' % (nodecnt, obj, solcnt))
            '''
            print(x)
            x_flat = []
            for i in range(len(x)):
                if x[i] == 1.0:
                    j = i % 20 + 1
                    x_flat.append(j)
            print(x_flat)
            '''
            intermediate_.append((runt, x, objbnd, objgap, solcnt, nodecnt, objbst))
            #intermediate_.append((runt, x, nodecnt, objgap, solcnt, objbnd))
            #total_energy = calc_energy_full(n, wind_speeds, prob_l, x, ss_coef)
            #print(total_energy)

        except GurobiError as e:
            print('Error code ' + str(e.errno) + ": " + str(e))
    
    #else:
    #    print('Non solution')

def gurobi_MILP(n, wind_speeds, prob_l, m, aggregated_coef, timelim):
    mod = Model("mip1")
    N = range(n)

    X = mod.addVars(N, vtype=GRB.BINARY, name='x')

    try:

        # Create a new model

        # singles = [i for i in N]

        # pairs = [(i,j) for i in N for j in N]

        Y = mod.addVars(N, N, vtype=GRB.BINARY, name='y')

        obj = LinExpr();
        obj += quicksum((0.33 * (np.dot(prob_l, pow(wind_speeds, 3))) * X[i] for i in
                         N))  # - energy_coef[(i,j)] * Y[i,j]) for i in N for j in N)
        # obj+= quicksum(0.33* prob_l [l] * pow(wind_speeds[l],3) * X[i] for i in N for l in range(len(wind_speeds)))
        obj -= quicksum(aggregated_coef[(i, j)] * Y[i, j] for i in N for j in N)
        mod.setObjective(obj, GRB.MAXIMIZE);

        # cSum = LinExpr();
        # cSum= quicksum(X[i] for i in singles)

        mod.addConstr((X.sum('*') == m), "c0")

        for i in range(n):
            # cSum += X[i]
            for j in range(n):
                if (i, j) in Epsilon:
                    mod.addConstr(X[i] + X[j] <= 1, name='cEps' + str(i) + str(j))
                mod.addConstr(X[i] + X[j] - Y[i, j] <= 1, name='cY' + str(i) + str(j))
        # mod.addConstr((cSum <= m), name = "cSum")

        # mod.addConstr((cSum == 0), name = "cSum")
        mod.Params.timeLimit = timelim


        # Optimize model
        #mod._vars = mod.getVars()

        mod.optimize(call_back)
        print(sum([X[i].x for i in N]))

        print('Obj: %g' % mod.objVal)

        # Print number of solutions stored
        # nSolutions = mod.SolCount
        # print('Number of solutions found: ' + str(nSolutions))
        #
        # # Print objective values of solutions
        # for e in range(nSolutions):
        #     mod.setParam(GRB.Param.SolutionNumber, e)
        #     print('Intermediate solution: ')
        #
        #     print(X.Xn)
        #     #for e in N:
        #     #    if X[e].Xn > .9:
        #     print('%g ' % mod.PoolObjVal, end='')
        #     if e % 15 == 14:
        #         print('')
        # print('')


    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))

    except AttributeError:
        print('Encountered an attribute error')

    #mod.setParam(GRB.Param.SolutionNumber, 1)
    x_vals = [X[i].x for i in N]
    total_energy = calc_energy_full(n, wind_speeds, prob_l, x_vals, ss_coef)
    print(total_energy)
    return total_energy
    #return X


def qubo_fy(wind_speeds, prob_l, n, aggregated_coef, m):
    # scaling if needed due to embedding considerations
    scaling = 50  # 0
    arr_x = []
    exp_x = {}
    for i in range(n):
        arr_x.append(Binary('x' + str(i)))
        cur_coef = 0.33 * np.dot(prob_l, pow(wind_speeds, 3))

        # scaling:
        # cur_coef = scaling*(0.33*(pow(u_inf,3)))/max_energy

        # objective value for setting a turbine in position i
        exp_x[i] = cur_coef * arr_x[i]
        # print(cur_coef)

    for i in range(n):

        if i % 50 == 0:
            print(str(i) + ' out of ' + str(n))

        for j in range(n):
            # only basic velocity componenet
            if i != j:
                # no scaling:
                cur_coef = aggregated_coef[i, j]
                # scaling:
                # cur_coef = scaling*energy_coef[i,j]/max_energy
                # the influence of (i,j) coefficient in the objective
                exp_x[i] -= cur_coef * arr_x[i] * arr_x[j]
    exp = 0
    for i in range(n):
        exp += exp_x[i]
    # setting a penalty term
    lam_ = 10500
    # computing sum of decision variables
    sum_x = sum([arr_x[i] for i in range(n)])

    # constraint on the number of turbines
    exp2 = (sum_x - m) ** 2

    exp3 = 0
    for i in range(n):
        for j in range(n):
            # must not locate turbines too closely
            if (i, j) in Epsilon:
                exp3 += arr_x[i] * arr_x[j]

    # total objective function (minimizing)
    exp_total = -exp + lam_ * (exp2) + lam_ * (exp3)

    # qubo encoding into xQx-c
    qub_exp = exp_total.compile().to_qubo()
    return qub_exp
# In[11]:

def check_matrix_sym(aggregated_coef):
    for i in range(n):
        for j in range(n):
            if i!=j:
                if aggregated_coef[i,j]!=aggregated_coef[j,i]:
                    print('assymetry detected')


def qubo_fy_eldan(wind_speeds, prob_l, n, aggregated_coef, m):
    matrix = np.zeros((n,n))
    # scaling if needed due to embedding considerations
    scaling = 50  # 0
    arr_x = []
    X = Array.create("x", n, 'BINARY')

    for i in range(n):
        matrix[i,i] = 0.33 * np.dot(prob_l, pow(wind_speeds, 3))


    for i in range(n):
        for j in range(i, n):
            matrix[i,j] -= aggregated_coef[i,j]
            matrix[j,i] -= matrix[i,j]


    lam_ = 10500
    # computing sum of decision variables

    for (i, j) in Epsilon:
        matrix[i, j] -= lam_



    expr = lam_ * (sum(X) - m) ** 2
    # constraint on the number of turbines
    constr_coeffs, constr_constant = expr.compile().to_qubo()
    for i in range(n):
        for j in range(i, n):
            val = max(constr_coeffs.get((X[i].label, X[j].label), 0), constr_coeffs.get((X[j].label, X[i].label), 0))
            matrix[i, j] -= val

    return matrix, constr_constant





### Gurobi version here.
# Benchmarking against MILP (Gurobi)
def gurobi_QUBO(n, wind_speeds, prob_l, m, aggregated_coef, timelim):
    N = range(n)
    mod = Model("mip1")
    X = mod.addVars(N, vtype=GRB.BINARY, name='x')

    try:
        lam_2 = 11000
        lam_1 = 1000
        # Create a new model

        # singles = [i for i in N]

        # pairs = [(i,j) for i in N for j in N]


        obj = QuadExpr();
        obj += quicksum((0.33 * (np.dot(prob_l, pow(wind_speeds, 3))) * X[i] for i in
                         N))  # - energy_coef[(i,j)] * Y[i,j]) for i in N for j in N)
        # obj+= quicksum(0.33* prob_l [l] * pow(wind_speeds[l],3) * X[i] for i in N for l in range(len(wind_speeds)))
        obj -= quicksum(aggregated_coef[(i, j)] * X[i]*X[j] for i in N for j in N)
        obj -= lam_1 * quicksum(X[i] * X[j] for i in N for j in N if (i, j) in Epsilon)
        obj -= lam_2 * (X.sum('*') - m) * (X.sum('*') - m)
        mod.setObjective(obj, GRB.MAXIMIZE);

        # cSum = LinExpr();
        # cSum= quicksum(X[i] for i in singles)

        # mod.addConstr((X.sum('*') == m), "c0")

        # for i in range(n):
        # cSum += X[i]
        #    for j in range(n):
        #        if (i,j) in Epsilon:
        #            mod.addConstr(X[i]+X[j]<=1, name = 'cEps'+str(i)+str(j))
        #        mod.addConstr(X[i]+X[j] - Y[i,j] <=1  , name = 'cY'+str(i)+str(j))
        # mod.addConstr((cSum <= m), name = "cSum")

        # mod.addConstr((cSum == 0), name = "cSum")
        mod.Params.timeLimit = timelim
        #setParam('OutputFlag', 0)
        #setParam('Heuristics', 0)
        # Optimize model
        mod.optimize(call_back)
        print(sum([X[i].x for i in N]))

        print('Obj: %g' % mod.objVal)

    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))

    except AttributeError:
        print('Encountered an attribute error')

    x_vals = [X[i].x for i in N]
    total_energy = calc_energy_full(n, wind_speeds, prob_l, x_vals, ss_coef)
    print(total_energy)
    return total_energy
# In[12]:
def gurobi_QILP(n, wind_speeds, prob_l, m, aggregated_coef, timelim):
    N = range(n)

    # Create a new model
    mod = Model("mip1")

    # singles = [i for i in N]

    # pairs = [(i,j) for i in N for j in N]

    X = mod.addVars(N, vtype=GRB.BINARY, name='x')
    try:


        obj = QuadExpr();
        obj += quicksum((0.33 * (np.dot(prob_l, pow(wind_speeds, 3))) * X[i] for i in
                         N))  # - energy_coef[(i,j)] * Y[i,j]) for i in N for j in N)
        obj -= quicksum(aggregated_coef[(i, j)] * X[i]*X[j] for i in N for j in N)
        mod.addConstr((X.sum('*') == m), "c0")

        for i in range(n):
            # cSum += X[i]
            for j in range(n):
                if (i, j) in Epsilon:
                    mod.addConstr(X[i] + X[j] <= 1, name='cEps' + str(i) + str(j))
        mod.setObjective(obj, GRB.MAXIMIZE);

        # cSum = LinExpr();
        # cSum= quicksum(X[i] for i in singles)

        # mod.addConstr((X.sum('*') == m), "c0")

        # for i in range(n):
        # cSum += X[i]
        #    for j in range(n):
        #        if (i,j) in Epsilon:
        #            mod.addConstr(X[i]+X[j]<=1, name = 'cEps'+str(i)+str(j))
        #        mod.addConstr(X[i]+X[j] - Y[i,j] <=1  , name = 'cY'+str(i)+str(j))
        # mod.addConstr((cSum <= m), name = "cSum")

        # mod.addConstr((cSum == 0), name = "cSum")
        mod.Params.timeLimit = timelim
        # Optimize model
        mod.optimize(call_back)
        print(sum([X[i].x for i in N]))

        print('Obj: %g' % mod.objVal)

    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))

    except AttributeError:
        print('Encountered an attribute error')
    x_vals = [X[i].x for i in N]
    total_energy = calc_energy_full(n, wind_speeds, prob_l, x_vals, ss_coef)
    print(total_energy)
    return total_energy
# In[12]:

def calc_SS_energy(n, wind_speeds, prob_l, X, ss_coef, U):
    total_energy = 0
    Energy = {}
    N = range(n)
    for i in N:
        for l in range(len(wind_speeds)):
            # only wake!!!
            Energy[(i, l)] = 0.33 * prob_l[l] * X[i].x
            temp = 0
            for j in N:
                if j in U[(i, l)]:
                    temp += X[j].x * ss_coef[i,j,l] #pow(1 - velocity[i, j, l] / wind_speeds[l], 2)
            Energy[(i, l)] *= pow(wind_speeds[l] * (1 - pow(temp, 0.5)), 3)
    total_energy = sum([Energy[(i, l)] for i in N for l in range(len(wind_speeds))])
    print(total_energy)
    return total_energy

# In[34]:



def plt_intermediate(intermediate_, timelim, da_sol, bottom, top):
    #initiating for the 3 models:
    X = []
    Y = []

    temp_X = []
    temp_Y = []
    last_t = intermediate_[0][0]

    for i in intermediate_:

        if last_t>i[0]:
            X.append(temp_X)
            Y.append(temp_Y)
            temp_X = [i[0]]
            temp_Y = [i[1]]
        else:
            temp_X.append(i[0])
            temp_Y.append(i[1])
        last_t = i[0]
    X.append(temp_X)
    Y.append(temp_Y)
    X.append(list(range(0,int(timelim)+1)))
    Y.append([0]*3 + [da_sol]*(int(timelim)-2))
    plt.style.use('seaborn-darkgrid')
    # create a color palette


    palette = plt.get_cmap('Set1')

    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    plt.rc('font', size=MEDIUM_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    num = 0
    # multiple line plot
    lab = ['MILP', 'QILP', 'QUBO', 'DA']
    for i in range(len(X)):

        plt.plot(X[i], Y[i], marker='', color=palette(num),linestyle = '-',
                 linewidth=2, label=lab[i])
        num+=1
    # Add legend
    plt.legend(loc=5, ncol=2)


    plt.ylim((bottom, top))
    # Add titles
    plt.title("Performance over time", loc='left', fontsize=12, fontweight=2, color='orange')
    plt.xlabel("Time (Seconds)")
    plt.ylabel("Power (kW)")
    plt.show()


    

#@Users: this is what you want to run with arguments from command prompt.
n_wind = int(sys.argv[1])
axis_ = int(sys.argv[2])
m = int(sys.argv[3])
timelim = int(sys.argv[4]) # set this to 36000 as 10 hours when analyzing problem difficulty

intermediate_ = []
n = pow(axis_, 2)  # cells
if axis_ == 10:
    step = 200  # cell size in meters
elif axis_ == 20:
    step = 100
else:
    step = 100

angle = 360 / n_wind

print('Calculating locations...')
location_x, location_y = calc_location(axis_, step)
print('Setting wind regime...')

wind_speeds, prob_l = set_wind_regime(n_wind)

print('Computing coefficients...')

aggregated_coef, ss_coef, Epsilon, U = \
    calc_coefficients(location_x, location_y, wind_speeds, prob_l, n)
print('Coefficients computed')
energies = []
config = []
method_names = []

prev_size = 0

df = pd.DataFrame(columns=['Runtime', 'Energy', 'Obj', 'Bound', 'Gap', 'SolCount', 'Node'])
df['Runtime'] = 0
df['Energy'] = 0
df['Obj'] = 0
df['Bound'] = 0
df['Gap'] = 1
df['SolCount'] = 0
df['Node'] = 0


print('Running Gurobi MILP...')
energies.append(gurobi_MILP(n, wind_speeds, prob_l, m, aggregated_coef, timelim))

# Uncomment the following two lines to run Gurobi QP
# print('Running Gurobi QP...')
# energies.append(gurobi_QILP(n, wind_speeds, prob_l, m, aggregated_coef, timelim))

# Uncomment the following two lines to run Gurobi QUBO
# print('Running Gurobi QUBO...')
# energies.append(gurobi_QUBO(n, wind_speeds, prob_l, m, aggregated_coef, timelim))

for j, i in enumerate(intermediate_):
    energy_calc = calc_energy_full(n, wind_speeds, prob_l, i[1], ss_coef)
    df = df.append({'Runtime': i[0], 'Energy': energy_calc, 'Obj': i[6], 'Bound': i[2], 'Gap': i[3], 'SolCount': i[4], 'Node': i[5]}, ignore_index=True)

filename = 'MP_' + str(n_wind) + '_' + str(axis_) + '_' + str(m) + '_' + str(timelim) + '.csv'
df.to_csv(filename, sep = ",", index = False)
