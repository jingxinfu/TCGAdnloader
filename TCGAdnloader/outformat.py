import os 

def storeData(df, parental_dir, sub_folder, cancer,index=True):
        store_folder = '/'.join([parental_dir, sub_folder])
        if not os.path.exists(store_folder):
                os.makedirs(store_folder)

        df.to_csv('/'.join([store_folder, cancer]), sep='\t', index=index)
