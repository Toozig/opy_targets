import pandas as pd
import requests
import json
from multiprocessing import Pool

DISEASE_GENE_QUERY = """
        query DiseaseAssociationsQuery($efoId: String!, $index: Int!, $size: Int!, $filter: String, $sortBy: String!, $aggregationFilters: [AggregationFilter!]) {
          disease(efoId: $efoId) {
            id
            associatedTargets(page: {index: $index, size: $size}, orderByScore: $sortBy, BFilter: $filter, aggregationFilters: $aggregationFilters) {
              count
              rows {
                target {
                  id
                  approvedSymbol
                  approvedName
                  __typename
                }
                score
                datatypeScores {
                  componentId: id
                  score
                  __typename
                }
                __typename
              }
              __typename
            }
            __typename
          }
        }
    """

SNP_DATA_QUERY = """
          query OpenTargetsGeneticsQuery(
          $ensemblId: String!
          $efoId: String!
          $size: Int!
        ) {
          disease(efoId: $efoId) {
            id
            evidences(
              ensemblIds: [$ensemblId]
              enableIndirect: true
              size: $size
              datasourceIds: ["ot_genetics_portal"]
            ) {
              rows {
                id
                disease {
                  id
                  name
                }
                diseaseFromSource
                studyId
                studySampleSize
                variantId
                variantRsId
                literature
                publicationYear
                publicationFirstAuthor
                pValueExponent
                pValueMantissa
                oddsRatio
                oddsRatioConfidenceIntervalLower
                oddsRatioConfidenceIntervalUpper
                beta
                betaConfidenceIntervalLower
                betaConfidenceIntervalUpper
                variantFunctionalConsequence {
                  id
                  label
                }
                resourceScore
                projectId
              }
            }
          }
        }

        """

VARIANT_QUERY= """
    query variantInfo($variantId: String!) {
        variantInfo(variantId: $variantId) {
            id
            rsId
            chromosome
            position
            refAllele
            altAllele
            nearestGene{
                symbol
            }
            nearestGeneDistance
            nearestCodingGeneDistance
            mostSevereConsequence
            caddRaw
            caddPhred
            gnomadAFR
            gnomadAMR
            # ... (other gnomad population fields)
        }
    }
    """
    

VARIANT_ID_QUERY = """
    query searchRsId($rsId: String!) {
        search(queryString: $rsId) {
            variants {
                id
            }
        }
    }
    """


def flatten_json(y):
    """
    Take a json and faltten it. 
    taken from the web
    """
    out = {}
  
    def flatten(x, name =''):
          
        # If the Nested key-value 
        # pair is of dict type
        if type(x) is dict:
              
            for a in x:
                flatten(x[a], name + a + '_')
                  
        # If the Nested key-value
        # pair is of list type
        elif type(x) is list:
              
            i = 0
              
            for a in x:                
                flatten(a, name + str(i) + '_')
                i += 1
        else:
            out[name[:-1]] = x
  
    flatten(y)
    return out


def ask_api(query_string, variables):
    """
    Recive query adn return a JSON string
    """

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})
#     print(r.status_code)

    # Transform API response into JSON 
    return json.loads(r.text)
    
def get_gene_data(target):
    """
    creates a sieries from a json of a SNP 
    """
    g_data = flatten_json(target['target'])
    for i in target['datatypeScores']:
        if i['componentId'] == 'genetic_association':
            g_data.update({'genetic_association_score': i['score']})
            continue
    return pd.Series(g_data).drop('__typename')

def get_disease_targets(efoId: str,sort_by: str='genetic_association', size: int=50, threshold=0.0035)-> pd.DataFrame:
    """
    efoId - ID of disease. example: ENSG00000196208 (endometriosis)
    sort_by - default: 'genetic_association'. Can be 'drugs' or any other option from opentarget
    size - amoun t of gene retrived by the query
    threshold - genetic association minimal score
    returns pd.DataFrame of the first -size- gene associated to the disease specified by the -efoId- sorted by -sortby-
    
    """

    v_dict = {
        "efoId": efoId, #disease ID
        "index": 0,
        "size": size,
        "sortBy": sort_by,
        "filter": "",
        "aggregationFilters": []
    }
    targets = ask_api(DISEASE_GENE_QUERY, v_dict)['data']['disease']['associatedTargets']['rows']
    result = pd.concat([get_gene_data(i) for i in targets],axis=1).T
    result = result[result.genetic_association_score >= threshold]
    return result

def get_snp_data(ens_id: str, asso_disease_id: str, size: int=10)-> pd.DataFrame:
    """
    ens_id - Ensembl_ID of gene. example: ENSG00000196208 (endometriosis)
    asso_disease_id - EFO ID of disease. example: EFO_0001065 (endometriosis)
    size - maximal amount of SNP retrived by the query
    returns DataFrame of all SNPs related to the gene and the given disease
    """
    var = {'ensemblId': ens_id,
     'efoId' :asso_disease_id, 
        'size' : size}

    api_response_as_json = ask_api(SNP_DATA_QUERY, var)
    snps = api_response_as_json['data']['disease']['evidences']['rows']
    
    if not len(snps):
        return pd.DataFrame()
    snps = [pd.Series({'variantRsId':snp['variantRsId'],
                'variantId': snp['variantId'],
               'studyId': snp['studyId'],
               'studySampleSize': snp['studySampleSize'],
              'publicationFirstAuthor': snp['publicationFirstAuthor'],
             'label': snp['variantFunctionalConsequence']['label']}) for snp in snps]
    df = pd.concat(snps, axis=1).T
    df['chr'] = df['chr'] = pd.DataFrame(df['variantId'].str.split('_').tolist())[0]
    df['location']= pd.DataFrame(df['variantId'].str.split('_').tolist())[1]
    df['disease_ID'] = asso_disease_id
    df['gene_related'] = ens_id
    return df
     
def get_SNP_df(efoId: str, n_gene:int =25, sort_by: str='genetic_association', size: int=50)-> pd.DataFrame:
    """
    efoId - EFO_ID of disease. example: EFO_0001065 (endometriosis)
    n_gene - default: 25. maximal number of genes to retrieve thier SNP associate with the disease 
    sort_by - default: 'genetic_association'. Can be 'drugs' or any other option from opentarget
    size - amoun t of gene retrived by the query
    returns DataFrame of all the the snps related to the -efoId- disease, takes SNPs only the top -n_gene-
    related according to the -sortBy- parameter 
    
    """
    disease_df = get_disease_targets(efoId, sort_by=sort_by, size=size)
    snp_result  = pd.DataFrame()
    for snp_id in disease_df['id'][:n_gene]:
        snp_result = pd.concat([snp_result, get_snp_data(snp_id,efoId)])
    return snp_result[~snp_result['variantRsId'].duplicated()].reset_index(drop=True)


# Function to search for a variant by its rsID using the Open Targets Genetics API
def search_variant_by_rsId(rsId):
    """
    Searches for a variant by its rsID using the Open Targets Genetics API.

    Args:
        rsId (str): The rsID of the variant to search for.

    Returns:
        str: The ID of the found variant or an empty string if not found.
    """
    url = 'https://api.genetics.opentargets.org/graphql'
    

    
    variables = {'rsId': rsId}
    
    search_result = requests.post(url, json={'query': VARIANT_ID_QUERY, 'variables': variables}).json()
    variants = search_result.get('data', {}).get('search', {}).get('variants', [])
    
    if variants:
        variant_id = variants[0]['id']
        return variant_id
    else:
        return ""

# Function to retrieve detailed information about a variant by its variant ID
def get_variant_info(variant_id):
    """
    Retrieves detailed information about a variant by its variant ID using the Open Targets Genetics API.

    Args:
        variant_id (str): The ID of the variant to retrieve information for.

    Returns:
        dict: A dictionary containing information about the variant.
    """
    url = 'https://api.genetics.opentargets.org/graphql'
    

    if variant_id == "":
        return ""
    
    variables = {'variantId': variant_id}
    
    variant_result = requests.post(url, json={'query': VARIANT_QUERY, 'variables': variables}).json()
    
    return variant_result['data']['variantInfo']

# Function to convert variant JSON data to a Pandas Series
def var_json_to_series(variant_json):
    """
    Converts variant JSON data to a Pandas Series.

    Args:
        variant_json (dict): A dictionary containing variant information.

    Returns:
        pd.Series: A Pandas Series containing variant information.
    """
    if variant_json == "":
        return pd.Series()
    
    variant = pd.Series(flatten_json(variant_json))  # Assuming flatten_json is defined elsewhere
    gnomAD_cols = variant.index[variant.index.str.contains('gnomad')]
    AF_popmax = variant[gnomAD_cols].max()
    variant = variant.drop(gnomAD_cols)
    variant['gnomAD2_AF_popmax'] = AF_popmax
    variant = variant.rename({'chromosome': 'CHROM', 'position': 'POS', 'refAllele': 'REF', 'altAllele': 'ALT'})
    return variant

# Function to get detailed information for a single variant by its rsID
def get_single_variant_info(rsid: str):
    """
    Gets detailed information for a single variant by its rsID.

    Args:
        rsid (str): The rsID of the variant.

    Returns:
        pd.Series: A Pandas Series containing detailed variant information.
    """
    var_id = search_variant_by_rsId(rsid)
    if var_id:
        var_json = get_variant_info(var_id)
        return var_json_to_series(var_json)

    
# Function to get detailed information for multiple variants
def get_variants_info(variant_series: pd.Series):
    """
    Gets detailed information for multiple variants.

    Args:
        variant_series (pd.Series): A Pandas Series containing variant rsIDs.

    Returns:
        pd.DataFrame: A Pandas DataFrame containing detailed information for the variants.
    """
    with Pool(processes=8) as pool:  # Adjust the 'processes' parameter based on your CPU cores
        results = pool.map(get_single_variant_info, variant_series)
    
    return pd.concat(results, axis=1).T.sort_values(['CHROM', 'POS', 'REF', 'ALT']).set_index(['CHROM', 'POS', 'REF', 'ALT'])
