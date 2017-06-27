# deseq output

def _generate_diff_expression_csv(self, result_directory, alpha_cutoff, fold_change_cutoff,
                                  condition_string):
    """
    _generate_diff_expression_csv: get different expression matrix with DESeq2
    """
    result_files = os.listdir(result_directory)
    if 'gene_count_matrix.csv' not in result_files:
        raise ValueError('Missing gene_count_matrix.csv, available files: {}'.format(
            result_files))

    rcmd_list = ['Rscript', os.path.join(os.path.dirname(__file__), 'run_DESeq.R')]
    rcmd_list.extend(['--result_directory', result_directory])
    rcmd_list.extend(['--alpha_cutoff', alpha_cutoff])
    rcmd_list.extend(['--fold_change_cutoff', fold_change_cutoff])
    rcmd_list.extend(['--condition_string', condition_string])

    rcmd_str = " ".join(str(x) for x in rcmd_list)

    self._run_command(rcmd_str)


def _generate_expression_matrix_file(self, result_directory):
    """
    _generate_expression_matrix_file: generate expression matrix file
    """

    expression_matrix_csv_file = os.path.join(result_directory, 'sig_genes_results.csv')
    expression_matrix_tsv_file = os.path.join(result_directory, 'sig_genes_results.tsv')
    with open(expression_matrix_csv_file, 'rb') as source:
        rdr = csv.reader(source)
        with open(expression_matrix_tsv_file, 'wb') as result:
            for r in rdr:
                result.write('\t'.join(r) + '\n')

    return expression_matrix_tsv_file


def _save_expression_matrix(self, result_directory, filtered_expression_matrix_name,
                            workspace_name):
    """
    _save_expression_matrix: save ExpressionMatrix object to workspace
    """
    log('start saving ExpressionMatrix object')

    expression_matrix_file = self._generate_expression_matrix_file(result_directory)

    tsv_file_to_matrix_params = {'input_file_path': expression_matrix_file,
                                 'genome_ref': self.expression_set_data.get('genome_id'),
                                 'data_type': 'log2_level',
                                 'data_scale': '1.0',
                                 'output_ws_name': workspace_name,
                                 'output_obj_name': filtered_expression_matrix_name}

    matrix_ref = self.fv.tsv_file_to_matrix(tsv_file_to_matrix_params)['output_matrix_ref']

    return matrix_ref
