import requests
import json
import pandas as pd
import re
import os

#os.chdir('rds/rds-cew54-basis/Projects/myositis-IMD/code')

# NB: Based on the sample script provided by Open Targets Genetics here: https://genyetics-docs.opentargets.org/data-access/graphql-api#available-endpoints

index_variants_and_studies_query = """
    query annotateLeadSnp($inputVariantId: String!) {
        indexVariantsAndStudiesForTagVariant(variantId: $inputVariantId) {
            associations {
            indexVariant {
            nearestGene {
                id
                symbol
            }
            mostSevereConsequence
            id
            }
            study {
                source
                traitReported
                pmid
                pubDate
                pubTitle
                pubAuthor
                hasSumstats
                nInitial
                nReplication
                nCases
                traitCategory
                numAssocLoci
            }
            pval
            beta
        }
  }
}
"""

variant_query = """
    query annotateLeadSnp($inputVariantId: String!) {
        variantInfo(variantId: $inputVariantId) {
            rsId
            nearestGene{
                id
                symbol
                bioType
            }
            nearestGeneDistance
            mostSevereConsequence
            }
        }
"""

genes_for_variant_query = """
query annotateLeadSnp($inputVariantId: String!){
    genesForVariant(variantId: $inputVariantId) {
        gene {
        id
        symbol
        bioType
        description
        }
        functionalPredictions {
        typeId
        sourceId
        aggregatedScore
        tissues {
            tissue {
            id
            }
            maxEffectLabel
            maxEffectScore
        }
        }
        distances {
        sourceId
        aggregatedScore
        tissues {
            tissue {
            id
            name
            }        
            distance 
        }
        }
        intervals {
        typeId
        sourceId
        tissues {
            score
        }
        }
        qtls {
        sourceId
        typeId
        aggregatedScore
                tissues {
            tissue {
            id
            }
                beta
                pval
                }
        }
    }
    }
"""

base_url = "https://api.genetics.opentargets.org/graphql"

daf = pd.read_csv('../data/snp.to.map.tsv', delimiter = '\t')

# daf.rename(columns = {'SNPID': 'SNP'}, inplace = True)

result_dict = {}

for index, row in daf.iterrows():

    # if re.match('^rs', row.SNP):
    #     variables = {'inputVariantId': row.SNP} # This script now takes rsids
    if re.match('\d+:\d+:\w+:\w+', row.SNP):
        variables = {"inputVariantId": row.SNP.replace(':', '_')}
    else:
        variables = {"inputVariantId": '_'.join([str(x) for x in [row.CHR, row.BP, row.REF, row.ALT]])}

    r = requests.post(base_url, json={"query": variant_query, "variables": variables})

    variant_response_data = json.loads(r.text)['data']['variantInfo']

    r = requests.post(base_url, json={"query": index_variants_and_studies_query, "variables": variables})

    index_variants_and_studies_response_data = json.loads(r.text)['data']['indexVariantsAndStudiesForTagVariant']

    r = requests.post(base_url, json={"query": genes_for_variant_query, "variables": variables})

    genes_for_variant_response_data = json.loads(r.text)['data']['genesForVariant']

    result_dict[row.SNP] = {'variantInfo' : variant_response_data, 'indexVariantsAndStudiesForTagVariant': index_variants_and_studies_response_data, 'genesForVariant': genes_for_variant_response_data}

with open('../data/mapped.genes_v2.tsv', 'w') as f:
    json.dump(result_dict, f)

d = []

for k,v in result_dict.items():
    d.append(
        {
            'SNPID': k,
            'rsID': v['variantInfo']['rsId'],
            'nearestGene': v['variantInfo']['nearestGene']['symbol'] if v['variantInfo']['nearestGene'] else ''
        }
        )


meta_daf = pd.DataFrame(d)

meta_daf = meta_daf.merge(right = daf, how = 'right', left_on = 'SNPID', right_on = 'SNP')

meta_daf.to_csv('../data/mapped.genes_v2.tsv', sep = '\t', index = False)
