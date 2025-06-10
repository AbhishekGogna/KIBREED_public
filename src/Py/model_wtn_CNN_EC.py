from .libs import *
from .func import CNN_net_fixed, CNN_net_flex, dense_post_concat

# Model 3 ---------------------------------------------------------------------------------------------------------------------
class tuner(kt.HyperModel):
    def __init__(self, marker_n, ec_n, tune = False, batch_norm = False):
        self.marker_n = marker_n # markers
        self.ec_n = ec_n
        self.tune = tune
        self.batch_norm = batch_norm
        print(f'markers = {self.marker_n} and g_a_s\'s = {self.ec_n}')

    def build(self, hp):
        # generate network for ec data
        net_g_s = CNN_net_fixed(hp, n_in = self.ec_n, 
                               model_name = "g_s", 
                               kernel_size = 2)
        # generate network for g_a data
        net_g_a = CNN_net_flex(hp, n_in = self.marker_n, 
                               model_name = "g_a", 
                               base_layer = 7, max_layer = 3, tuning = self.tune)
        #net_g_d = CNN_g_data(hp, self.marker_n, "g_d")
        X_list = [net_g_a]
        X_list.append(net_g_s) # should follow the order of inputs provided
        
        # split list for input and output
        input_list = [x[0] for x in X_list]
        output_list = [x[1] for x in X_list]
        
        # concat inputs
        model_concat = concatenate(output_list, name = "concat_in")
        
        if self.batch_norm:
            # do normalization
            model_concat = BatchNormalization(name = "norm_in")(model_concat) # i don't know enoough about this to be confident in using it
        
        model_compiled = dense_post_concat(hp, input_list, model_concat, tuning = not self.tune)
        
        return model_compiled
