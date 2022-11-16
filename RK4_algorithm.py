import numpy as np



#Fourth order Runge-Kutta algorithm optimized for sparse matices
class RK4_algorithm:
    def __init__(self, start, end, iterations, commConst, operator, initial_state, eqsNum, sparse_matrix_indexes, input_labels, output_labels, norm_prec = 0.0):
        self.start = start
        self.end = end
        self.iterations = iterations
        self.commConst = commConst
        self.operator = operator
        self.c = initial_state
        self.eqsNum = eqsNum
        self.dt = (end-start)/(iterations-1)
        self.sparse_matrix_indexes = sparse_matrix_indexes
        self.input_label = input_labels
        self.output_labels = output_labels

    # The execution of the fourth order Runge-Kutta algorithm
    def calculate(self):
        print('calc')
        times = np.linspace(self.start, self.end,self.iterations) # időlépések

        dt_p2 = self.dt/2

        num_of_steps = -1


        for time in times:
            num_of_steps+=1
            prev_step_num = num_of_steps-1
            if num_of_steps == 0:
                continue

            if num_of_steps%(np.round(self.iterations*0.01,0)) == 0:
                print("Progress: " + str(num_of_steps/self.iterations*100))
                #Checking the norm of the present coefficients of the atom's states
                coeffs = list(map(lambda x: abs(x)**2, self.c[:,-1]))
                curr_norm = sum(coeffs)
                print(curr_norm)


            slope_a = []
            slope_b = []
            slope_c = []
            slope_d = []


            last_element = -1 # index of the last element


            # The self.spM.trip object helps the algorithm to travel only to the non-zero elements of the matice
            # This part calculates the slopes
            for row in self.sparse_matrix_indexes.trip:
                slope_a.append(self.create_part_slope_a( time,row, prev_step_num))

            for row in self.sparse_matrix_indexes.trip:
                slope_b.append(self.create_part_slope_b( time, row, prev_step_num, slope_a, dt_p2))

            for row in self.sparse_matrix_indexes.trip:
                slope_c.append(self.create_part_slope_c( time, row, prev_step_num, slope_b, dt_p2))

            for row in self.sparse_matrix_indexes.trip:
                slope_d.append(self.create_part_slope_d( time, row, prev_step_num, slope_c))


            #After calculating the slopes it aggregates them
            curr_res = []
            for state_num in range(0, self.eqsNum):
                curr_res.append([0])
                curr_res[state_num][0] = self.c[state_num][last_element] + self.dt / 6 * (slope_a[state_num] + 2 * slope_b[state_num] + 2 * slope_c[state_num] + slope_d[state_num])

            #The current state is inserted
            self.c = np.append(self.c,curr_res,axis=1)

        self.result = { self.input_label : times }
        for row, i in zip(self.c, range(0,len(self.c))):
            self.result[self.output_labels[i]] = row
        return self.result



    # These functions calculates the slopes
    def create_part_slope_a(self ,time,spMRow, prev_step_num):
        part_slope_a =0
        for element in spMRow:
            i = element[0]
            j = element[1]
            part_slope_a += self.commConst * self.operator(time,i,j) * self.c.item((j,prev_step_num))        
        return part_slope_a

    def create_part_slope_b(self,  time, spMRow, prev_step_num, slope_a, dt_p2):
        part_slope_b = 0
        for element in spMRow:
            i = element[0]
            j = element[1]            
            part_slope_b += self.commConst * self.operator(time + dt_p2,i,j) * (self.c.item((j,prev_step_num)) + dt_p2 * slope_a[j])
        return part_slope_b

    def create_part_slope_c(self, time, spMRow, prev_step_num, slope_b, dt_p2):
        part_slope_c = 0
        for element in spMRow:
            i = element[0]
            j = element[1]            
            part_slope_c += self.commConst * self.operator(time + dt_p2,i,j) * (self.c.item((j,prev_step_num)) + dt_p2 * slope_b[j])
        return part_slope_c
    
    def create_part_slope_d(self,  time, spMRow, prev_step_num, slope_c):
        part_slope_d = 0
        for element in spMRow:
            i = element[0]
            j = element[1]  
            part_slope_d += self.commConst * self.operator(time + self.dt,i,j) * (self.c.item((j,prev_step_num)) + self.dt * slope_c[j])
        return part_slope_d