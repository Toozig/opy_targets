import pandas as pd
import requests
import json

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
