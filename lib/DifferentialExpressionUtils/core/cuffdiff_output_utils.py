# cuffdiff output


def generate_and_upload_expression_matrix(logger, scratch, rscripts, scriptfile, shock_url, hs_url, token, cuffdiff_dir,
                                          ws_url, workspace):
    TSV_to_FeatureValue = "trns_transform_TSV_Exspression_to_KBaseFeatureValues_ExpressionMatrix"
    returnVal = False
    if exists(cuffdiff_dir) == False:
        logger.info("Cuffdiff directory does not exist")
        return False

    # generate expression matrix
    # input = Rscript fpkmgenematrix.R cuffdiff_dir outpath(os.join)

    outmatrix = join(scratch, scriptfile) + ".matrix"
    # outmatrixparse =  join (scratch, scriptfile) + ".matrix.parse.txt"
    outmatrix2 = scriptfile + ".matrix.txt"
    outjson = scriptfile + ".matrix.parse.txt.json"

    computescript = join(rscripts, scriptfile)
    if (exists(computescript) == False):
        logger.info("Rscript does not exist")
        return False

    # Generate command to be executed
    ropts = ["Rscript", computescript]

    ropts.append("--cuffdiff_dir")
    ropts.append(cuffdiff_dir)

    ropts.append("--out")
    ropts.append(outmatrix)

    roptstr = " ".join(str(x) for x in ropts)

    # Run Rscript to generate Expression matrix
    openedprocess = subprocess.Popen(roptstr, shell=True, stdout=subprocess.PIPE)
    openedprocess.wait()

    if openedprocess.returncode != 0:
        logger.info("R script did not return normally, return code -"
                    + str(openedprocess.returncode))
        return False

    matrix_parse = parse_expression_matrix_separate_comma(outmatrix)
    # convert expression matrix TSV to json

    cmd_expression_json = [TSV_to_FeatureValue,
                           '--workspace_service_url', ws_url,
                           '--workspace_name', workspace,
                           '--object_name', matrix_parse,
                           '--working_directory', scratch,
                           '--input_directory', scratch,
                           '--output_file_name', outjson]

    cmd_expression_json = ["perl", '/kb/module/lib/kb_cummerbund/get_exp_matrix.pl',
                           '/kb/module/work/outmatrix.parse.txt', '/kb/module/work/out.json']
    logger.info(" ".join(cmd_expression_json))
    tool_process = subprocess.Popen(" ".join(cmd_expression_json), stderr=subprocess.PIPE, shell=True)
    stdout, stderr = tool_process.communicate()

    if stdout is not None and len(stdout) > 0:
        logger.info(stdout)
    if stderr is not None and len(stderr) > 0:
        logger.info(stderr)

    if tool_process.returncode != 0:
        return False

    return outjson


def filter_expression_matrix(fparams, system_params):
    cuffdiff_dir = fparams['cuffdiff_dir']
    scratch = system_params['scratch']
    selected_condition_option = fparams['pairs']
    sample1 = fparams['sample1']
    sample2 = fparams['sample2']
    q_value_cutoff = abs(float(fparams['q_value_cutoff']))
    log2_fold_change_cutoff = abs(float(fparams['log2_fold_change_cutoff']))
    print q_value_cutoff
    print log2_fold_change_cutoff
    infile = fparams['infile']
    outfile = fparams['outfile']
    num_genes = 1000000000000000000000  # no upper bound
    try:
        num_genes = int(fparams['num_genes'])
    except:
        num_genes = 100
    outf = f.filter_expresssion_matrix_option(scratch, infile, outfile, sample1, sample2, num_genes, q_value_cutoff,
                                              log2_fold_change_cutoff)
    if (os.path.isfile(outf)):
        return outf
    else:
        return False


def create_expression_matrix(self, ctx, expressionMatrixParams):
    """
    :param expressionMatrixParams: instance of type
       "expressionMatrixParams" -> structure: parameter "workspace_name"
       of type "workspace_name" (workspace name of the object), parameter
       "ws_cuffdiff_id" of type "ws_cuffdiff_id" (@id ws
       KBaseRNASeq.RNASeqCuffdiffdifferentialExpression), parameter
       "ws_expression_matrix_id" of type "ws_expression_matrix_id" (@id
       ws KBaseFeatureValues.ExpressionMatrix), parameter
       "include_replicates" of type "bool" (indicates true or false
       values, false <= 0, true >=1)
    :returns: instance of type "ws_expression_matrix_id" (@id ws
       KBaseFeatureValues.ExpressionMatrix)
    """
    # ctx is the context object
    # return variables are: returnVal
    #BEGIN create_expression_matrix

    params    = expressionMatrixParams
    returnVal = params['ws_expression_matrix_id']
    #Set up workspace client
    user_token = ctx['token']
    workspace = params['workspace_name']
    ws_client  = Workspace(url=self.__WS_URL, token=user_token)

    #Read the input cuffdiff workspace object json file and get filehandle for cuffdiff tar file
    s_res = ws_client.get_objects([{
        'name' : params['ws_cuffdiff_id'],
        'workspace' : params['workspace_name']
        }])

    # Check if workspace has data
    if len(s_res) == 0:
        self.__LOGGER.info("Workspace did not return any objects")
        return returnVal

    cuffdiff_dir = join (self.__SCRATCH , "cuffdiffData/cuffdiff")
    cuffdiff_dir = script_util2.extract_cuffdiff_data(self.__LOGGER, self.__SHOCK_URL, self.__SCRATCH, s_res, user_token)
    self.__LOGGER.info("Cuffdiff folder = " + cuffdiff_dir)

    if (cuffdiff_dir is False):
        return returnVal

    # Run R script to get fpkmgenematrix.R

    # Prepare output object.
    outjson = False;
    #outjson = "repfpkmgenematrix.R.matrix.txt.json";

    if params['include_replicates'] ==0:
     scriptfile = "fpkmgenematrix.R"
     outjson = script_util2.generate_and_upload_expression_matrix(self.__LOGGER, self.__SCRATCH,
                self.__RSCRIPTS, scriptfile, self.__SHOCK_URL, self.__HS_URL, user_token,
                cuffdiff_dir, self.__WS_URL,workspace)


    else:
     scriptfile = "repfpkmgenematrix.R"
     outjson = script_util2.generate_and_upload_expression_matrix(self.__LOGGER, self.__SCRATCH,
                self.__RSCRIPTS, scriptfile, self.__SHOCK_URL, self.__HS_URL, user_token,
                cuffdiff_dir, self.__WS_URL,workspace)

    if outjson is False:
        self.__LOGGER.info("Creation of expression matrix failed")
        return returnVal
    with open("{0}/{1}".format(self.__SCRATCH , outjson),'r') as et:
              eo = json.load(et)
    eo['type']='untransformed'
    genome_ref = s_res[0]['data']['genome_id']
    eo['genome_ref'] = genome_ref

    self.__LOGGER.info(workspace + self.__SCRATCH + outjson + params['ws_expression_matrix_id'])
    ws_client.save_objects({'workspace' : workspace,
        'objects' : [{ 'type' : 'KBaseFeatureValues.ExpressionMatrix',
                       'data' : eo,
                       'name' : params['ws_expression_matrix_id']
                    }]})

    #END create_expression_matrix

    # At some point might do deeper type checking...
    if not isinstance(returnVal, basestring):
        raise ValueError('Method create_expression_matrix return value ' +
                         'returnVal is not type basestring as required.')
    # return the results
    return [returnVal]

