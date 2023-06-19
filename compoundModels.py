import logging, sys

#####--------------------COMPOUND MODEL START-------------------------------#####

class ISP_Compound():
    def __init__(self, goodputRatio = 1, **segments):
        self.segments = segments
        self.goodputRatio = goodputRatio
        self.min_throughput = -1
        self.bottleNeck = ""
        
    def __str__(self):
        map = "Client X"
        for segmentName, segmentObject in self.segments.items():
            map += "-----"+segmentName+": "+ segmentObject.__str__() + "-----X"
        map += " Server"
        return map

    def compute_throughput(self):
        for segmentName, segmentObject in self.segments.items():
            segmentObject.avg_throughput()
            logging.debug(f"{segmentName} has an average throughput of {segmentObject.ssThroughput}")
            if (segmentObject.ssThroughput < self.min_throughput) or (self.min_throughput==-1):
                self.min_throughput = segmentObject.ssThroughput
                self.bottleNeck = segmentName
        return self.min_throughput
    
    def goodput(self):
        #Here I assume that the bottleneck is the satellite link
        if self.min_throughput == -1:
            self.compute_throughput()
        
        return self.min_throughput*self.goodputRatio
    
    def time_to_transfer(self,filesize=100, goodput = True):
        # Still need to add qpep delay
        if goodput:
            return filesize/self.goodput()
        return filesize/self.min_throughput