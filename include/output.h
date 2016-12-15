#ifndef OUTPUT_H
#define OUTPUT_H

#include "simulator.h"
#include "pedigree.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <htslib/vcf.h>

namespace sim {
namespace io {


unsigned int allele_orders[4][3] = {{1, 2, 3}, {2, 3, 0}, {3, 0, 1}, {0, 1, 2}};


//TODO: Much of the functionality of VCFOutput and TadOutput can be moved into a parent class.
class VCFOutput {
 protected:
  htsFile *vcf_handle_;
  bcf_hdr_t *vcf_hdr_;
  bcf1_t *vcf_rec_;

  std::string c_contig_id_;
  Base c_ref_;
  int32_t c_pos_;

  
 public:
  VCFOutput(std::string &filename) {
    VCFOutput(filename.c_str());
  }

  VCFOutput(const char *filename) {
    vcf_handle_ = hts_open(filename, "w");
    vcf_hdr_ = bcf_hdr_init("w");
    vcf_rec_ = bcf_init();
    bcf_hdr_append(vcf_hdr_, "##fileformat=VCFv4.2");
    bcf_hdr_append(vcf_hdr_, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  }

  void addContig(const char *contig, int length) {
    std::string contig_header = std::string("##contig=<ID=") + contig +",length=" + std::to_string(length) + ">";
    bcf_hdr_append(vcf_hdr_, contig_header.c_str());
  }

  void addSample(const char *sample) {
    bcf_hdr_add_sample(vcf_hdr_, sample);
  }
  
  void writeHeader() {
    bcf_hdr_add_sample(vcf_hdr_, nullptr);
    bcf_hdr_write(vcf_handle_, vcf_hdr_);
  }

  void newSite(std::string &contig, int32_t pos, Base ref) {
    c_contig_id_ = contig;
    c_pos_ =  pos;
    c_ref_ = ref;
  }

  
  void addGenotypes(std::vector<Genotype> &gts) {
    if(c_ref_ == N) {
      vcf_rec_->pos = c_pos_;
      bcf_update_alleles_str(vcf_hdr_, vcf_rec_, "N");
      std::vector<int32_t> genotypes(2*bcf_hdr_nsamples(vcf_hdr_), bcf_gt_unphased(0));
      bcf_update_genotypes(vcf_hdr_, vcf_rec_, &genotypes[0], genotypes.size());	
      return;
    }
    
    std::string alleles = std::string(base2str(c_ref_)); // string of alleles
    std::array<int, 4> allele_map = {-1,-1,-1,-1};
    allele_map[c_ref_] = 0;
    int allele_map_last = 1;
    
    std::vector<int32_t> genotypes(2*gts.size());
    for(int a = 0; a < gts.size(); ++a) {
      Genotype gt = gts[a];
      Base b1 = split_table[gt][0];
      int32_t b1_index = allele_map[b1];
      if(b1_index == -1) {
	b1_index = allele_map[b1] = allele_map_last++;
	alleles += std::string(",") + base2str(b1); 
      }

      Base b2 = split_table[gt][1];
      int32_t b2_index = allele_map[b2];
      if(b2_index == -1) {
	b2_index = allele_map[b2] = allele_map_last++;
	alleles += std::string(",") + base2str(b2); 
      }

      if(b1_index <= b2_index) {
	genotypes[a*2] = bcf_gt_unphased(b1_index);
	genotypes[a*2+1] = bcf_gt_unphased(b2_index);
      }
      else {
	// NOT sure if the standard requires the reference goes first for unphases genotypes, code may be unnnecessary.
	genotypes[a*2] = bcf_gt_unphased(b2_index);
	genotypes[a*2+1] = bcf_gt_unphased(b1_index);
      }      
    }

    vcf_rec_->pos = c_pos_;
    bcf_update_alleles_str(vcf_hdr_, vcf_rec_, alleles.c_str());
    bcf_update_genotypes(vcf_hdr_, vcf_rec_, &genotypes[0], genotypes.size());					

  }

  void writeSite() {
    bcf_write(vcf_handle_, vcf_hdr_, vcf_rec_);
  }
  
  void close() {
    bcf_hdr_destroy(vcf_hdr_);
    bcf_close(vcf_handle_);
  }
  
};

 
class TadOutput {
 protected:
  std::vector<std::pair<std::string, int>> contigs_;
  std::vector<std::string> libraries_;

  Base c_ref_;
  std::string c_contig_;
  size_t c_pos_;
  std::vector<reads_list> c_depths_;

  std::ofstream outfile_;
  
  
 public:
  TadOutput(std::string &filename) {
    outfile_.open(filename.c_str());
  }

  TadOutput(const char *filename) {
    outfile_.open(filename);
  }

  void addContig(std::string &contig_name, int len) {
    contigs_.emplace_back(contig_name, len);
    
  }
  
  void addLibrary(std::string &lb_name) {
    libraries_.emplace_back(lb_name);    
    
  }

  void writeHeader() {
    outfile_ << "@ID\tFF:TAD\tVN:0.1\tSO:coordinate" <<  std::endl;
    for(std::pair<std::string, int> &p : contigs_) {
      outfile_ << "@SQ\tSN:" << p.first << "\tLN:" << std::to_string(p.second) << std::endl;
    }
    
    for(std::string &m : libraries_) {
      outfile_ << "@AD\tID:" << m << "\tSM:" << m << std::endl;
    }
    
  }

  void newSite(std::string &contig_name, size_t pos, Base ref) {
    c_contig_ = contig_name;
    c_pos_ = pos;
    c_ref_ = ref;
  }

  void addCalls(std::vector<reads_list> &depths) {
    c_depths_ = depths;
  }
  
 
  void writeSite() {
    // Get references
    if(c_ref_ == N) {
      // Don't display N's, maybe later.
      //outfile_ << c_contig_ << "\t" << c_pos_ << "\tN" << std::endl;
      return;
    }
    
    std::string alleles = base2str(c_ref_);
    std::vector<std::string> depths_str(libraries_.size());
    for(int a = 0; a < c_depths_.size(); ++a) {
      depths_str[a] = std::to_string(c_depths_[a][c_ref_]);
    }

    for(Base allele : allele_orders[c_ref_]) {
      int total = 0;
      for(int lb = 0; lb < c_depths_.size(); ++lb) {
	total += c_depths_[lb][allele];
      }
      if(total > 0) {
	alleles += base2str(allele);
	for(int lb = 0; lb < c_depths_.size(); ++lb) {
	  depths_str[lb] += "," + std::to_string(c_depths_[lb][allele]);
	}
      }      
    }
    
    outfile_ << c_contig_ << "\t" << c_pos_ << "\t" << alleles;
    for(std::string &dp : depths_str) {
      outfile_ << "\t" << dp;
    }
    outfile_ << std::endl;
      
  }

  void close() {
    outfile_.close();
  }
  
};


class CSVOutput {
 private:
  std::ofstream outfile_;

 public:
  CSVOutput(const char *filename) {
    outfile_.open(filename);
    outfile_ << "ID,father,mother,"
	     << "Ns,ref/homo,ref/het,-/hom,-/het," // Info on diff between founder and reference
	     << "Mendelian match,1 Mut (Germ),2 Mut (Germ)," // Number of different germline mutations
	     << "Somatic match,1 Mut (Som),2 Mut (Som)" // Number of different somatic mutations
	     << std::endl;

  }


  void addStats(std::vector<ped::Member> &family) {
    for(ped::Member &m : family) {
      outfile_ << m.id << ",";
      if(m.is_founder) {
	outfile_ << ".,.,"
		 << m.stats.n_unknowns << ","<< m.stats.n_hom_ref_match << "," << m.stats.n_het_ref_match << "," << m.stats.n_hom_snp << "," << m.stats.n_het_snp << ","
		 << ".,.,.,";
	  
      }
      else {
	outfile_ << get_dad(m).id << "," << get_mom(m).id << ","
		 << ".,.,.,.,.,"
		 << m.stats.n_mendelians << "," << m.stats.n_single_mut << "," << m.stats.n_double_mut << ",";
      }

      outfile_ <<  m.stats.n_som_matches << "," << m.stats.n_som_single_mut << "," << m.stats.n_som_double_mut << std::endl;
    }
  }
  
  void close() {
    outfile_.close();
  }

  
};
 
} // namespace io
} // namespace sim

#endif
