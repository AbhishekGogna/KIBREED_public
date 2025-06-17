from .libs import *
from .func import CNN_net_fixed, CNN_net_flex, dense_post_concat

# Model 3 ---------------------------------------------------------------------------------------------------------------------
class tuner(kt.HyperModel):
    def __init__(self, marker_n, ev_n, tune = False):
        self.marker_n = marker_n # markers
        self.ev_n = ev_n
        self.tune = tune
        print(f'markers = {self.marker_n} and env variables = {self.ev_n}')

    def build(self, hp):
        # generate network for ev data
        net_ev = CNN_net_fixed(hp, n_in = self.ev_n, 
                               model_name = "ev", 
                               kernel_size = 2)
        # generate network for g_a data
        net_g_a = CNN_net_flex(hp, n_in = self.marker_n, 
                               model_name = "g_a", 
                               base_layer = 7, max_layer = 3, tuning = self.tune)
        X_list = [net_g_a]
        X_list.append(net_ev) # should follow the order of inputs provided
        
        # split list for input and output
        input_list = [x[0] for x in X_list]
        output_list = [x[1] for x in X_list]
        
        # concat inputs
        model_concat = concatenate(output_list, name = "concat_in")
        
        model_compiled = dense_post_concat(hp, input_list, model_concat, tuning = not self.tune)
        
        return model_compiled
