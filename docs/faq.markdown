---
layout: page
title: FAQ
permalink: /faq/
---

## Frequently Asked Questions


#### **1. How many sequences can RECUR handle?**
 
   The RECUR algorithm is fast and light weight and can be used on tens of thousands of sequences. However, if a constraint tree is not provided, RECUR is limited by the scalability of IQ-TREE. 

#### **2. Can I use IQ-TREE3?**
   The short answer is, yes you can!
   
   RECUR requires IQ-TREE version 2.0 or higher to function correctly. The default IQ-TREE version bundled with RECUR is 2.4.0, which uses the binary name `iqtree2`. If you're using a different version of IQ-TREE (e.g., IQ-TREE3), ensure that the corresponding binary is accessible in your environment. You can specify the binary name using the `-iv` or `--iqtree-version` flag when running RECUR. For example, to use IQ-TREE3, add `-iv iqtree3` to your command.

   Different versions of IQ-TREE may have variations in their command-line interfaces. To accommodate these differences, RECUR provides a [config.json](https://github.com/OrthoFinder/RECUR/releases/download/v1.0.0/config.json) file that allows users to define the appropriate command structure for their specific IQ-TREE version.
   
   If a new version of IQ-TREE (e.g., IQ-TREE4) is released, you can adapt the config.json by copying an existing configuration block (such as for `iqtree2` or `iqtree3`) and renaming the method key to match the new binary name (iqtree4). Adjust the command parameters as needed to align with the new version's requirements.

   When running RECUR with a custom `config.json`, specify its path using the `--config` flag and indicate the IQ-TREE binary with the `-iv` flag:
   ```bash
   recur [options] -f <alignment_file> --outgroups <outgroup_species/file> -st <AA|CODON> -te <treefile> --config path/to/config.json -iv iqtree4
   ```
   This ensures that RECUR utilizes the correct command structure for the specified IQ-TREE version.​

   ⚠️ Caution: While minor adjustments to the command structure can be managed via the `config.json`, significant changes in IQ-TREE's interface may not be compatible with RECUR. In such cases, it's recommended to contact the RECUR developers for guidance on updating the package accordingly.

---

*Still have questions? Please head to our [dicussion section](https://github.com/OrthoFinder/RECUR/issues) on GitHub to check for the similiar questions and answers. If no solutions were found, please open an new issue.*