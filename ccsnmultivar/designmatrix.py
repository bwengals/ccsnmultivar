# design matrix class and methods
import re
import warnings
import numpy as np
import csv
import itertools


class DesignMatrix(object):
    def __init__(self,formula):
        self._formula = formula
        # parse formula
        self._formula_dict, self._inter_list = _parse_formula(formula)
        self._trimmed_columns = []

    def get_metadata(self):
        # return column names, row names, formula
        print "not done yet"

    def get_rownames(self):
        return self._row_names

    def get_matrix(self):
        # actually make the design matrix
        X_array, _ = self.run_encoder(self._param_dict, self.encoder)
        return X_array

    def get_formula(self):
        return self._formula

    def trim_columns(self, columns_to_trim):
        """
        remove column in design matrix
        """
        # TODO check if trimmed column is actually one of the columns
        if len(self._trimmed_columns) == 0:
            self._trimmed_columns.append(columns_to_trim)
        else:
            self._trimmed_columns.extend(columns_to_trim)
        self._trimmed_columns = self._trimmed_columns[0]
        self.encoder['trimmed_columns'] = self._trimmed_columns

    def undo_column_trims(self):
        """
        add back in the columns that were trimmed by calls to trim_columns()
        """
        self._trimmed_columns = []
        self.encoder['trimmed_columns'] = self._trimmed_columns

    def get_columnnames(self):
        # get all columns in full columns list that havent been trimmed
        columns = self._full_columns
        columns_to_trim = self._trimmed_columns
        [columns.remove(item) for item in columns_to_trim]
        return columns


    def make_param_dict_from_file(self,path_to_params):
        """
        make param dict from a file on disk
        """
        # then we were given a path to a parameter file
        param_list = list(csv.reader(open(path_to_params,"rb")))
        # delete empty elements (if any)
        param_file = [x for x in param_list if x != []]
        # make dict of [wavenames] = raw_params
        name_list = []
        param_list = []
        # get header names for each param (names of param_list columns)
        param_colnames = param_file[0][1:] # 0th element is "Name" or "Wavename"
        # start from 1.  (row 0 is the header)
        for i in np.arange(1, len(param_file)):
            name_list.append(param_file[i][0])
            param_list.append(param_file[i][1:])
        # remove ' ' blank spaces from param_list
        param_list = [[x.strip() for x in y] for y in param_list]
        param_dict = {}
        # i loops through param_colnames, j loops thru param values per wave
        for i in np.arange(0, len(param_colnames)):
            param_dict[param_colnames[i]] = []
            for j in np.arange(0,len(name_list)):
                param_dict[param_colnames[i]].append(param_list[j][i])
        # now we have param_dict, and name_list
        self._param_dict = param_dict
        self._row_names = name_list

    def __call__(self, *args):
        self.make_param_dict_from_file(args[0])

        # when design matrix is given parameters, make the encoder
        self.encoder    = self.make_encoder(self._formula_dict,self._inter_list,
                                            self._param_dict)
        # actually make the design matrix
        _, self._full_columns = self.run_encoder(self._param_dict, self.encoder)
        return self

    def make_encoder(self,formula_dict,inter_list,param_dict):
        """
        make the encoder function
        """
        X_dict = {}
        Xcol_dict = {}
        encoder_dict = {}
        # first, replace param_dict[key] = values, with param_dict[key] = dmatrix
        for key in formula_dict:
            encoding,arg = formula_dict[key]
            if 'Dev' in encoding:
                # make deviation encoded design matrix
                drop_name = arg
                # encode
                deviation_encoder,X_sub,colnames_sub = _dev_encode(param_dict,drop_name,key)
                # additionally, store in dictionary for use by interactions
                X_dict[key] = X_sub
                Xcol_dict[key] = colnames_sub
                # store dictionary of encoder functions to keep for prediction
                encoder_dict[key] = deviation_encoder
            elif 'Dum' in encoding:
                # make dummy variable encoding design mat
                ref_name = arg
                dummy_encoder,X_sub,colnames_sub = _dum_encode(param_dict,ref_name,key)
                # additionally, store in dictionary for use by interactions
                X_dict[key] = X_sub
                Xcol_dict[key] = colnames_sub
                # store dictionary of encoder functions to keep for prediction
                encoder_dict[key] = dummy_encoder
            elif 'Poly' in encoding:
                # make polynomial encoding design mat
                degree = arg
                polynomial_encoder,X_sub,colnames_sub = _poly_encode(param_dict,degree,key)
                # additionally, store in dictionary for use by interactions
                X_dict[key] = X_sub
                Xcol_dict[key] = colnames_sub
                # store dictionary of encoder functions to keep for prediction
                encoder_dict[key] = polynomial_encoder
            else:
                print encoding
                raise Exception("Encoding name error")
        # now compute interaction designmatrices
        for interaction in inter_list:
            if len(interaction) >= 3:
                raise Exception("Doesn't allow 4-way or higher interaction terms")
            elif len(interaction) == 3:

                param_name1 = interaction[0]
                param_name2 = interaction[1]
                param_name3 = interaction[2]
                col_names1 = Xcol_dict[param_name1]
                col_names2 = Xcol_dict[param_name2]
                col_names3 = Xcol_dict[param_name3]

                # make 3-way encoder function
                def threeway_encoder(param_name1,param_name2,param_name3, \
                                     col_names1, col_names2, col_names3, X_dict):
                    """
                    needs the three names of the parameters to be encoded, as well as
                    a dictionary containing the already encoded single parameter 
                    design matrices, keyed by name
                    """
                    X1 = X_dict[param_name1]
                    X2 = X_dict[param_name2]
                    X3 = X_dict[param_name3]

                    X_int = []
                    names_int = []
                    for i in np.arange(0,X1.shape[1]):
                        for j in np.arange(0,X2.shape[1]):
                            for k in np.arange(0,X3.shape[1]):
                                X_int.append(X1[:,i]*X2[:,j]*X3[:,k])
                                names_int.append(col_names1[i] + "*" + \
                                                 col_names2[j] + "*" + col_names3[k])
                    # make X_int from lists to np array
                    X_int = np.array(X_int).T
                    return X_int, names_int
                encoder_dict['threeway'] = threeway_encoder

            elif len(interaction) == 2:
                # there are two interaction terms (A*B)

                param_name1 = interaction[0]
                param_name2 = interaction[1]
                col_names1  = Xcol_dict[param_name1]
                col_names2  = Xcol_dict[param_name2]

                # make twoway_encoder function
                def twoway_encoder(param_name1,param_name2, col_names1, col_names2, X_dict):
                    X1 = X_dict[param_name1]
                    X2 = X_dict[param_name2]

                    X_int = []
                    names_int = []
                    for i in np.arange(0,X1.shape[1]):
                        for j in np.arange(0,X2.shape[1]):
                            X_int.append(X1[:,i]*X2[:,j])
                            names_int.append(col_names1[i] + "*" + col_names2[j])
                    X_int = np.array(X_int).T
                    return X_int, names_int
                encoder_dict['twoway'] = twoway_encoder
            else:
                raise Exception("Error while evaluating meaning of interaction term")

        # make key in encoder to specify which columns are active
        encoder_dict['trimmed_columns'] = self._trimmed_columns
        return encoder_dict


    def run_encoder(self,param_dict, encoder_dict):
        """
        run the encoder on a supplied param_dict
        """
        X_dict = {}
        Xcol_dict = {}
        # put each column of X in Xbycol_dict
        Xbycol_dict = {}
        for key in encoder_dict:
            if (key != 'twoway') and (key != 'threeway') and (key != 'trimmed_columns'):
                encoder = encoder_dict[key]
                param_values = param_dict[key]
                Xsub,names = encoder(key,param_values)
                X_dict[key] = Xsub
                Xcol_dict[key] = names
                for i in np.arange(0,len(names)):
                    Xbycol_dict[names[i]] = Xsub[:,i]

        # now do interactions
        inter_list = self._inter_list
        for interaction in inter_list:
            if 'twoway' in encoder_dict.keys():
                encoder = encoder_dict['twoway']
                param_name1 = interaction[0]
                param_name2 = interaction[1]
                col_names1 = Xcol_dict[param_name1]
                col_names2 = Xcol_dict[param_name2]
                X_int, names = encoder(param_name1,param_name2, \
                                       col_names1,col_names2, X_dict)

                # put columns into Xbycol_dict
                for i in np.arange(0,len(names)):
                    Xbycol_dict[names[i]] = X_int[:,i]

            if 'threeway' in encoder_dict.keys():
                encoder = encoder_dict['threeway']
                param_name1 = interaction[0]
                param_name2 = interaction[1]
                param_name3 = interaction[2]
                col_names1 = Xcol_dict[param_name1]
                col_names2 = Xcol_dict[param_name2]
                col_names3 = Xcol_dict[param_name3]
                X_int, names = encoder(param_name1,param_name2,param_name3, \
                                       col_names1,col_names2,col_names3, X_dict)

                # put columns into Xbycol_dict
                for i in np.arange(0,len(names)):
                    Xbycol_dict[names[i]] = X_int[:,i]

        # remove columns that were trimmed (if any)
        trimmed_columns = encoder_dict['trimmed_columns']
        full_columns = Xbycol_dict.keys()
        used_columns = [x for x in full_columns if x not in trimmed_columns]
        # make design matrix array
        X = []
        for name in used_columns:
            X.append(Xbycol_dict[name])
        # always add intercept column last
        X.insert(0,np.ones(np.shape(X[0])))
        used_columns.insert(0,'Intercept')
        # final design matrix
        X = np.vstack(X).T
        return X, used_columns




def _parse_formula(formula):
    """
    Parse formula into a dictionary
    formula_dict[variable_name] = [encoding, dropped_name]
    Parse interactions into a list
    inter_list = [[A,B], [A,C], [A,B,C]]

    formula = "A + beta + A*beta | Dev(A,drop=1), Poly(beta,degree=3)"
    """
    #TODO: DEAL WITH BSPLINE HAS MULTIPLE ARGUEMENTS (NOT JUST ONE)

    # break formula apart from encoding instructions
    formula,instr = formula.replace(' ','').split('|')
    # break formula apart at + sign
    formula_terms = formula.split('+') 
    # examine the instructions term, first split by ),
    instr = instr.split('),')
    # elements in the instr list will match the number of non interaction 
    #      elements in formula
    # go through each formula term, making a dictionary whose key is variable name
    formula_dict = {}
    encoding     = []
    other_arg    = []
    inter_list   = []
    for term in formula_terms:
        if "*" in term:
            # then this is an interaction term, make 'list of lists'
            term = term.split('*')
            inter_list.append(term)
        else:
            # then this is not an interaction term, make blank dictionary
            formula_dict[term] = [encoding, other_arg]
    # loop through instructions, parse each term
    for term in instr:
        # remove punctuation in term
        term = re.sub('[()]','',term)
        # check for each encoding type
        if "Dev" in term:
            # remove first three letters (Dev)
            term = term[3:]
            # split on comma
            var_name,arg = term.split(',')
            # split on equals sign in arg, take part after
            drop_name = arg.split('=')[1]
            # put var_name and drop_name into proper key in formula_dict
            formula_dict[var_name] = ["Dev",drop_name]
        elif "Dum" in term:
            # remove first three letters (Dum)
            term = term[3:]
            # split on comma
            var_name,arg = term.split(',')
            # split on equals sign in arg, take part after
            ref_name = arg.split('=')[1]
            # put var_name and drop_name into proper key in formula_dict
            formula_dict[var_name] = ["Dum",ref_name]
        elif "Poly" in term:
            # remove first four letters (Poly)
            term = term[4:]
            # split on comma
            var_name,arg = term.split(',')
            # split on equals sign in arg, take part after
            degree = arg.split('=')[1]
            # put var_name and drop_name into proper key in formula_dict
            formula_dict[var_name] = ["Poly",degree]
        else:
            raise Exception("Unknown Encoding")
    return formula_dict,inter_list





def _dev_encode(param_dict,drop_name,param_name):
    """
    _dev_encode takes the parameter dictionary in, as well as the name of parameter
    to drop from design matrix, and the string name of the parameter

    it returns the encoded design matrix, list of column name strings, and
    and encoder function to encode new values of param properly
    """

    # param_dict[col_name] = np.array(wave_values)
    values = param_dict[param_name]
    # treat as string, map each unique value to an integer
    unique_values = []
    [unique_values.append(i) for i in values if not unique_values.count(i)]
    # warning if number of unique values is large compared to number of waveforms
    if len(unique_values) > 0.5*len(values):
        warning.warn("Many unique values in parameter %s for a Deviation encoding" & param_name)
    # make deviation encoding instruction matrix
    I = np.eye(len(unique_values))
    # remove column of I for drop_name
    drop_idx = unique_values.index(drop_name)
    I = np.delete(I,drop_idx,1)
    # make row of I for drop_name index
    I[drop_idx,:] = -1.
    # make map_dict[unique_value[i]] = I_row
    map_dict = {}
    for i in np.arange(0,len(unique_values)):
        map_dict[unique_values[i]] = np.float64(I[i,:])

    # make parameter to designmatrix encoder function
    def deviation_encoder(param_name,new_values):
        # TODO: first check if new_values are outside the range of the originals
        #if set(new_values) != set(unique_values):
        #    raise ValueError("Cannot extrapolate outside original parameter range")
        X = np.empty((len(new_values),len(map_dict[new_values[0]]-1)))
        for i in np.arange(0,len(new_values)):
            X[i,:] = map_dict[new_values[i]]

        # make names out of values
        col_names = []
        for i in np.arange(0,len(unique_values)):
            col_names.append(param_name+":["+unique_values[i]+" - "+"mu"+"]")
        # remove column name for drop_name
        col_names.remove(param_name + ":[" + drop_name + " - " + "mu" + "]")

        return X, col_names

    # use deviation_encoder to encode values
    X, col_names = deviation_encoder(param_name,values)
    return deviation_encoder, X, col_names


def _dum_encode(param_dict,drop_name,param_name):
    # param_dict[col_name] = np.array(wave_values)
    values = param_dict[param_name]
    # get unique values in order of appearance
    unique_values = []
    [unique_values.append(i) for i in values if not unique_values.count(i)]
    # warning if number of unique values is large compared to number of waveforms
    if len(unique_values) > 0.5*len(values):
        warning.warn("Many unique values in parameter %s for a Deviation encoding" & param_name)
    # make deviation encoding instruction matrix
    I = np.eye(len(unique_values))
    # remove column of I for drop_name
    drop_idx = unique_values.index(drop_name)
    I = np.delete(I,drop_idx,1)

    # make map_dict[unique_value[i]] = I_row
    map_dict = {}
    for i in np.arange(0,len(unique_values)):
        map_dict[unique_values[i]] = np.float64(I[i,:])

    # use map_dict to encode design matrix
    X = np.empty((len(values),len(map_dict[unique_values[0]]-1)))
    for i in np.arange(0,len(values)):
        X[i,:] = map_dict[values[i]]

    # make parameter to designmatrix encoder function
    def dummy_encoder(param_name, new_values):
        # TODO: first check if new_values are outside the range of the originals
        #if set(new_values) != set(unique_values):
        #    raise ValueError("Cannot extrapolate outside original parameter range")    
        X = np.empty((len(new_values),len(map_dict[new_values[0]]-1)))
        for i in np.arange(0,len(new_values)):
            X[i,:] = map_dict[new_values[i]]
        # make column name strings
        col_names = []
        for i in np.arange(0,len(unique_values)):
            col_names.append(param_name+":["+unique_values[i]+" - " + drop_name + "]")
        # remove column name for drop_name
        col_names.remove(param_name + ":[" + drop_name + " - " + drop_name + "]")
        return X, col_names

    # use dummy_encoder to encode
    X, col_names = dummy_encoder(param_name,values)
    return dummy_encoder, X, col_names



def _poly_encode(param_dict,degree,param_name):
    x = np.array(map(float,param_dict[param_name]))
    degree = int(degree)
    if len(np.unique(x)) < 5:
        warning.warn("Not many unique values of parameter %s for Poly encoding" & param_name)
    n = degree
    min_x = min(x)
    max_x = max(x)

    # make encoder function
    def polynomial_encoder(param_name, new_values):  
        # generate chebyshev polynomials
        x = new_values
        x = np.array(map(float, x))
        # check if new_values are outside the range of the originals
        if (max(x) > max_x) or (min(x) < min_x):
            raise ValueError("Cannot extrapolate outside original parameter range")          
        m = len(x)
        # generate the z variable as a mapping of your x data into [-1,1]
        z = ((x - min_x)-(max_x - x))/(max_x - min_x)
        X = np.empty([m,n+1])
        X[:,0] = np.ones([m,])
        if n >= 1:
            X[:,1] = z;
        if n >= 2:
            for k in np.arange(2, n+1):
                X[:,k] = 2.*z*X[:,k-1] - X[:,k-2]
        # remove first column of A (never keep intercept)
        X = X[:,1:]
        # generate names for each column
        col_names = []
        for i in np.arange(0,degree):
            col_names.append(param_name + "^" + str(i+1))
        return X, col_names

    X, col_names = polynomial_encoder(param_name,x)
    return polynomial_encoder, X, col_names

def _bs_encode(param_dict,bsplineparameters, param_name):
    print "not implemented yet"
