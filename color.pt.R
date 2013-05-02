color.pt<-function(p.value,fc,min.p,min.fc)
{
	
	n<-length(p.value)
	color<-rep("lightgrey",n)
	for (i in 1:n)
	{
	if (p.value[i]>=min.p) color[i]="lightgrey"
	else
		{
		if(fc[i]>min.fc) color[i]="red"
		if(fc[i]<(-min.fc)) color[i]="green3"
		}
	}
	color
}
