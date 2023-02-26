import numpy as np
import pandas as pd
class GnomAD_df:
        
        
    
    def __init__(self,path, peak_file=None, remove_unkown=True, remove_phased_gt=True, only_peak_variants=True):
        """
        This class warps the result of gnomAD and interval  mATAC/hATAC analysis.
        path - path to the parquet dataframe file 
        peak_file - bed file which contains the following columns - CHROM,FROM,TO,INTERVAL_ID, ~BUT w/o header~
        remove_unknown - removed unkwon genotypes (e.g "./.")' 
        remove_phased_gt - "replaced phased genotype (e.g "0|1" -> "0/1")"
        only_peak_variants - removing variants outside of peak interval
        
        """
        self.__original_df = self.__open_table(path)
        self.f_phasing = False
        self.__filter_samples = []
        self.__filters_description = []
        self.__filter_functions = []
        self.__cur_df = self.__original_df.copy()
        self.__applied_filter_index = 0
        if only_peak_variants:
            self.remove_non_peak_variants()
        if remove_unkown:
            self.remove_unkown()
        if remove_phased_gt:
            self.remove_phasing()
        self.__samples = [i.split(':')[0] for i in self.__original_df.columns if i.endswith('GT')]
        if peak_file != None:
            self.__peak_df = self.__get_peak_df(peak_file)
    
    
    def __get_peak_df(self, path):
        peak_df = pd.read_csv(path, sep='\t', header=None)
        peak_df.columns = ['CHROM', 'FROM','TO', 'INTERVAL_ID']
        return peak_df
    
    def load_peak_file(peak_file_path):
        """
        loads a peak file into the object
        """
        self.__peak_df = self.__get_peak_df(peak_file)
    
    def get_peak_df(self):
        """
        returns the intervals dataframe 
        """
        return self.__peak_df
    
    def __open_table(self, path):
        df = pd.read_parquet(path).fillna(value=np.nan)
        interval = df.INTERVAL_ID.copy()
        df = df.drop(columns='INTERVAL_ID')
        df.insert(0,'INTERVAL_ID',interval)
        return df

    def __remove_non_peak_variants(self):
        """
        remove variants if they are not in a peak
        """
        self.__cur_df = self.__cur_df.replace('.', np.nan)
        self.__cur_df = self.__cur_df[~self.__cur_df.INTERVAL_ID.isna()]

    def remove_non_peak_variants(self):
        decs="removing variants outside of peak interval"
        
        if self.__check_if_filter_exists(decs):
            self.__filters_description.append("removing variants outside of peak interval")
            self.__filter_functions.append(self.__remove_non_peak_variants)
        return self

    def __remove_unknow(self):
        self.__cur_df = self.__cur_df.replace('./.', np.nan).replace('.|.', np.nan)

    def remove_unkown(self):
        """
        removed unkwon genotypes (e.g "./.")
        """
        decs='removed unkwon genotypes (e.g \"./.\")'
        
        if self.__check_if_filter_exists(decs):
            self.__filters_description.append('removed unkwon genotypes (e.g \"./.\")')
            self.__filter_functions.append(self.__remove_unknow)
        return self

    def remove_phasing(self):
        """
        "replaced phased genotype (e.g "0|1" -> "0/1")"
        """
        desc = "replaced phased genotype (e.g \"0|1\" -> \"0/1\")"
        if self.__check_if_filter_exists(desc):
            self.__filters_description.append(desc)
            self.__filter_functions.append(self.__remove_phasing)
        return self

    def __remove_phasing(self):

        """
        replace all "X|X" into "X/X"
        """
        if self.f_phasing: return
        for col in [i for i in self.__cur_df.columns if i.endswith('GT')]:
            for val in self.__cur_df[col].unique():
                if val == np.nan: continue
                if '|' in str(val):
                    self.__cur_df = self.__cur_df.replace(val, "%s/%s" % tuple(val.split('|')))
        self.f_phasing = True
    
    
    def __check_if_filter_exists(self,description):
        """
        Indicate if a filter already applied on the DF by checking if the description exist in the list.
        Return True if filter was not applied
        """
        if description in self.__filters_description:
            print(f"{description.split('(')[0]} already applied")
            return False
        return True
        
        
    
    def filter_dp(self, dp_t:float):
        """
        Fiter the variants (replace the GT with NaN) according to a given threshold (dp_t)
        """
        desc = f"Removing variants with reading depth {dp_t} or below"
        if self.__check_if_filter_exists(desc):
            self.__filters_description.append(desc)
            self.__filter_functions.append(lambda: self.__filter_dp(dp_t))
        return self

    def __filter_dp(self, dp_t:float):
        for sample in self.__samples:
            self.__cur_df.loc[self.__cur_df[sample + ':DP'] <= dp_t, sample + ':GT'] = np.nan

    def ___filter_AF_remove_unknown(self, af_t:float):
        self.__cur_df = self.__cur_df[(self.__cur_df.AF <= af_t) | (self.__cur_df.AF.isna())]

    def ___filter_AF(self,af_t:float):
        self.__cur_df = self.__cur_df[(self.__cur_df.AF <= af_t)]

    def filter_AF(self, af_t: float,remove_unkwon=False):
        """
        Filters the variants of gnomAD's allele frequency (AF) below a given threshold.
        If remove_unkown is True, remove variants with no record on gnomAD
        """
        description = f"Removing variants with allele frequency above {af_t}"
        description += " and variants with no record on gnomAD" if remove_unkwon else ""
        if self.__check_if_filter_exists(description):
            self.__filters_description.append(description)
            if not remove_unkwon:
                self.__filter_functions.append(lambda :self.___filter_AF_remove_unknown(af_t))
            else:
                self.__filter_functions.append(lambda :self.___filter_AF(af_t))
        return self

    def samples_filter(self,sample_list, appened=True):
        """
        Filter the variants dataframe to a given set of samples.
        If appened is False, remove all previous given samples of filteration

        """
        if appened:
            self.filter_samples += sample_list
        else:
            self.filter_samples = sample_list
        return self

    def __sample_filter(self):
        to_remove_samples = [i for i in self.__samples if i in self.__filter_samples]
        drop_cols = []
        for i in to_remove_samples:
            self.__cur_df = self.__cur_df.drop(columns=[f'{i}:GT', f'{i}:DP'])

    def print_filters(self):
        """
        prints the filters applied on the variant df
        """
        for f in range(len(self.__filters_description)):
            print(f"{f}. {self.__filters_description[f]}")

    def remove_filter(self, i):
        """
        remove the i'th filter from the table 
        To know the filters index, use "print_filter"
        """
        if i >= len(self.__filters_description):
            print("Filter does not exist")
            return self
        del self.__filters_description[i]
        del self.__filter_functions[i]
        if i <= self.__applied_filter_index:
            self.__cur_df = self.__original_df.copy()
            self.f_phasing = False
            self.__applied_filter_index = 0
        return self

    def reset_table(self, remove_phased_gt=True, remove_unkown=True, only_peak_variants=True):
        """
        returns the dataframe to it's initial state.
        removes all filters applied on the dataframe
        """
        self.__filter_samples = []
        self.__filters_description = []
        self.__filter_functions = []
        self.__applied_filter_index = 0
        self.__samples = [i.split(':')[0] for i in self.__original_df.columns if i.endswith('GT')]
        self.__cur_df = self.__original_df.copy()
        self.f_phasing = False
        if only_peak_variants:
            self.remove_non_peak_variants()
        if remove_phased_gt:
            self.remove_phasing()
        if remove_unkown:
            self.remove_unkown()
        return self

    def get_filters_decription(self):
        """
        returns list of filter description
        """
        return self.__filters_description
        
    def __apply_filters(self, verbos=True):
        max_index = len(self.__filters_description)
        # filter samples only if there are samples in the list
        if len(self.__filter_samples):
            if verbos:
                print("Filtering samples")
            self.__sample_filter()
        if verbos:
            for i in range(self.__applied_filter_index):
                print(self.__filters_description[i])
        for i in range(self.__applied_filter_index, max_index):
            if verbos:
                print(self.__filters_description[i])
            self.__filter_functions[i]()
        self.__applied_filter_index = max_index


    def get_table(self, verbos=True):
        """
        returns the variants df filtered by all given filters.
        Verbos = true -> prints the log of filteration
        """
        if verbos:
            print("applying filters")
        self.__apply_filters(verbos)
        if verbos:
            print('getting table')
        return self.__cur_df.copy()

    def __remove_reference(self):
        """
        removing 0/0 and 0/. or ./0 from the table (not saving the dataframe)
        """
        letters = ['0', '.']
        df = self.__cur_df.copy()
        for i in letters:
            for j in letters:
                df = df.replace(f"{i}|{j}", np.nan)
        return df

    def bool_variant_df(self):
        """
        returns a boolean table of samples VS variants.
        True if the variant exist in the sample.
        NOTICE- THIS APPLY REMOVE PHASING FILTER ON DF
        """
        self.__remove_phasing()
        self.__apply_filters()
        df = self.__remove_reference()
        gt = df[[i for i in df.columns if i.endswith('GT')]].notna()
        gt.columns = [i.replace(':GT','') for i in gt.columns ]
        return pd.concat([df[['AF','INTERVAL_ID']],gt], axis=1)
    
    
    
    def __add_interval_info(self,df,drop_all_zero=False):
        """
        assuming df index is the interval_id
        """
        peak_df_isin = self.__peak_df[self.__peak_df.INTERVAL_ID.isin(df.index)].set_index('INTERVAL_ID')
        result = pd.concat([peak_df_isin,df], axis=1)
        if not drop_all_zero:
            peak_df_isout = self.__peak_df[~self.__peak_df.INTERVAL_ID.isin(df.index)].set_index('INTERVAL_ID')
            result = pd.concat([result, peak_df_isout])
        return result
    
    
    
    def count_interval(self, boolean=False, drop_all_zero=False, per_sample=False):
        """
        returns a dataframe with the amount of variant each  in each peak.
        if boolean = True, reaturns a df of interval vs samples with Trued if sample ahve varint in that interval
        if drop_all_zero = True, drop all interval with zero variants 
        if per_sample = True returns a dataframe with the amount of variant per sample in each peak.
        """
        df =  self.bool_variant_df().drop(columns=['AF'])
        sum_df = df.groupby('INTERVAL_ID').sum()
        sum_df = sum_df if not boolean else sum_df != 0
        if not per_sample:
            sum_df = sum_df.sum(axis=1)
        return self.__add_interval_info(sum_df,  drop_all_zero)      
    
    def count_peaks_variants(self,numeric=False):
        """
        returns data frame with the amount of people who have variants in each interval
        if numeric is True - the total number of variants there are in each peak of all samples
        """
        bool_df = self.bool_variant_df().drop(columns=['AF'])
        melt_df = bool_df.melt(id_vars=['INTERVAL_ID'])
        melt_df = melt_df[melt_df.value]
        if numeric:
            melt_df = melt_df.groupby('INTERVAL_ID')
        else:
            melt_df = melt_df.drop_duplicates().groupby('INTERVAL_ID')
        result = melt_df.value.sum().sort_values()
        if numeric:
            result = result.rename(columns={0:'n_variants'})
        return 
    
    def sum_peaks(self):
        """
        sum the amount of variants in each peak 
        """
        bool_df = self.bool_variant_df()
        return bool_df.groupby('INTERVAL_ID').sum()
    
    def get_samples(self):
        return [i for i in self.__samples if i not in self.__filter_samples]
    
    def add_ken_to_main(self, to_add_path, sample_name=False, header=True, ken=False):
        """
        This functions add kens tsv data to the main df
        ken - the data is kens source
        """
    
        main_df = self.__original_df.copy()
        if not sample_name:
            sample_name = to_add_path.split('/')[-1].split('.')[0]
        if sample_name in self.__samples:
            print(f"{sample_name} already exist")
            return self
        to_add = pd.read_csv(to_add_path, sep='\t', header=None if not header else 0)
        col_names = [f'{sample_name}:GT',f'{sample_name}:DP']
        columns = ['CHROM','POS', 'REF','ALT','INTERVAL_ID','AF'] + col_names
        to_add.columns = columns
        if ken:
            to_add.ALT = to_add.ALT.str.replace(',<NON_REF>','')
        to_add = to_add.replace('.',np.nan).set_index(['CHROM', 'POS', 'ALT', 'REF']).fillna(value=np.nan)
        to_add[f'{sample_name}:DP'] = to_add[f'{sample_name}:DP'].astype(float)
        to_add.AF = to_add.AF.astype(float)
        exist = to_add[to_add.index.isin(main_df.index)]
        non_exist = to_add[~to_add.index.isin(main_df.index)]
#         return main_df,exist[col_names]
        main_df = pd.concat([main_df,exist[col_names]],axis = 1)

        main_df = pd.concat([main_df, non_exist.drop(columns=col_names)])
        main_df.loc[non_exist.index,col_names] = non_exist[col_names]
        self.__original_df = main_df
        self.reset_table()
        return self
    
    def count_sample_variants(self,verbos=True):
        """
        Counts the amount of variant each sample have.
        """
        bool_df = self.bool_variant_df(verbos=verbos)
        return bool_df.drop(columns=['AF','INTERVAL_ID']).sum()
    
    def save_df(self,save_format, path, original=False):
        """
        save the df in the format mentioned: csv/prq(parquet)
        """
        df = self.__original_df if original else self.get_table() 
        if save_format == 'csv':
            df.to_csv(f"{path}.{save_format}")
        elif  save_format == 'prq':
            df.to_parquet(f"{path}.{save_format}")
    
    def concat(gnomAD_df_list: list, reset_df=False, verbos=True):
        """
        Generate a new gnomAD_df from a list of gnomAD_df.
        The function takes the filtered DF from each object and creates a new one.
        all previous filters will not be transfered, hence- the basiv DF of the new df
        won't have the original's dfs data
        """
        df_list = []
        template = 'for samples %s the following filteres were applied:'
        for df in gnomAD_df_list:
            if verbos:
                samples = ",".join(df.get_samples())
                print(template % samples)
                for i in df.get_filters_decription():
                    print(i)
            df_list.append(df.get_table(verbos = False))
        new_df = pd.concat(df_list, axis=1)
        return GnomAD_df(new_df)
    
    

                      